#ifndef POSITION_H
#define POSITION_H

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map> 
#include <functional>
#include "AES.h"

#include <chrono>
#include <cmath>
#include <bitset>
#include <ctime>

#include <omp.h>

#include "NTL/ZZX.h"

#include "helib/FHE.h"
#include "helib/EncryptedArray.h"
#include "helib/replicate.h"

using namespace std;

extern NTL::ZZX long2Poly(long query);

class Node {
    public:
    map<string, int> childs;
    string keyword;
    // vector<NTL::ZZX> bitmap_plaintext;
    // vector<Ctxt *> bitmap_cipher;
    // Ctxt* cipher;
    // NTL::ZZX en_bitmap;
    int ciphertext_len;

    Node(int bitmap_size, int ptxt_num, const FHEPubKey& publicKey){
        ciphertext_len = 0;
        // bitmap_plaintext.resize(bitmap_size * ptxt_num, NTL::ZZX(0));
        // bitmap_cipher = new Ctxt(publicKey);
        
        // for (int i = 0; i < ptxt_num; i++) {
        //     bitmap_cipher.push_back(new Ctxt(publicKey));
        // }
        // cipher = new Ctxt(publicKey);
    }

    Node (const Node& p){
        childs = p.childs;
        keyword = p.keyword;
        // bitmap = p.bitmap;
        ciphertext_len = p.ciphertext_len;
        // bitmap_plaintext = p.bitmap_plaintext;
        // for (int i = 0; i < p.bitmap_cipher.size(); i++) {
        //     bitmap_cipher.push_back(new Ctxt(*p.bitmap_cipher[i]));
        // }
        // cipher = new Ctxt(*p.cipher);
        // bitmap_cipher = NULL;
        // cout << "2:" << bitmap_cipher << endl;
    }

    Node (Node&& p) {
        // cout << "3" << endl;
        childs = p.childs;
        keyword = p.keyword;
        // bitmap = p.bitmap;
        ciphertext_len = p.ciphertext_len;
        // bitmap_plaintext = p.bitmap_plaintext;
        // bitmap_cipher = p.bitmap_cipher;
        // for (int i = 0; i < p.bitmap_cipher.size(); i++)
        //     p.bitmap_cipher[i] = NULL;
        
        // cipher = p.cipher;
        // p.cipher = NULL;
    }

    Node& operator=(const Node& p){
        // cout << "4" << endl;
        childs = p.childs;
        keyword = p.keyword;
        // bitmap = p.bitmap;
        ciphertext_len = p.ciphertext_len;
        // bitmap_plaintext = p.bitmap_plaintext;
        // for (int i = 0; i < p.bitmap_cipher.size(); i++) {
        //     bitmap_cipher.push_back(new Ctxt(*p.bitmap_cipher[i]));
        // }
        // cipher = new Ctxt(*p.cipher);
        return *this;
    }

    Node& operator=(Node&& p){
        // cout << "5" << endl;
        if (this != &p){
            childs = p.childs;
            keyword = p.keyword;
            // bitmap = p.bitmap;
            ciphertext_len = p.ciphertext_len;
            // bitmap_plaintext = p.bitmap_plaintext;

            // for(int i = 0; i < this->bitmap_cipher.size(); i++)
            //     delete bitmap_cipher[i];

            // bitmap_cipher = p.bitmap_cipher;
            // for (int i = 0; i < p.bitmap_cipher.size(); i++)
            //     p.bitmap_cipher[i] = NULL;

            // cipher = p.cipher;
            // p.cipher = NULL;
        }
        return *this;
    }

    ~Node(){
        // cout << "~Node: " << bitmap_cipher << endl;
        // for(int i = 0; i < bitmap_cipher.size(); i++)
        //     delete bitmap_cipher[i];
        // delete cipher;
    }

    // int encryptBitmap(EncryptedArray &ea, const FHEPubKey &publicKey){
    //     ea.encrypt(*bitmap_cipher, publicKey, bitmap_plaintext);
    //     return 1;
    // }

    // int decryptBitmap(EncryptedArray &ea, FHESecKey &secretKey, vector<NTL::ZZX> &plaintextResult){
    //     ea.decrypt(*bitmap_cipher, secretKey, plaintextResult);
    //     return 1;
    // }
};


class PositionHeap {
    char *T;
    int n;
    // FHESecKey *secretKey_pointer;
    // EncryptedArray *ea_pointer;
    // FHEcontext *context_pointer;
    
public:
    vector<Node> nodes;
    vector<Ctxt> cipher;
    // void test_operator(){
    //     cout << "start test_operator()\n";
    //     cout << "context_pointer: " << context_pointer << endl;
    //     *(nodes[0].bitmap_cipher) += *(nodes[18].bitmap_cipher);
    //     cout << "end test_operator()\n";
    // }
    // T: "ab_ab#cd_cd#" 
    // keywords_list: {"ab", "cd"}
    PositionHeap(const char T[], 
                //  unsigned char *aes_key, 
                 vector<vector<string>> &records, 
                 int attribute_index,
                 FHESecKey *secretKey_pointer,
                 EncryptedArray *ea_pointer,
                 FHEcontext *context_pointer,
                 int dualposition_flag) {

        const FHEPubKey &publicKey = *secretKey_pointer;
        // unsigned char iv[AES_BLOCK_SIZE];
        int row_ind = records.size();
        int ptxt_num = ceil((float)row_ind / ea_pointer->size());
        int copy_flag;
        n = (int)strlen(T);
        int cipher_num = ceil((float) (n + 1) / ea_pointer->size());
        this->T = strdup(T);
        cout << "n:" << n << endl;

        nodes.resize(n+1, Node(ea_pointer->size(), ptxt_num, publicKey));
        cipher.resize(cipher_num, Ctxt(publicKey));

        // vector<NTL::ZZX> bitmap_plaintext(ea_pointer->size() * ptxt_num, NTL::ZZX(0));
        vector<vector<NTL::ZZX>> plain;
        for (int i = 0; i < cipher_num; i++) {
            vector<NTL::ZZX> temp_plain;
            temp_plain.resize(ea_pointer->size(), NTL::ZZX(0));
            plain.push_back(temp_plain);
        }
        // ea_pointer->encrypt(cipher[0], publicKey, plain);
        // cout << "after nodes.resize\n";
        int flag_count = 0;
        
        for (int i = n; --i >= 0; ) {

            const char *q = &T[i];
            if ((*q) == '_'){  // '_' links a keyword with its cope
                copy_flag = 0; //  
                continue;
            }else if ((*q) == '#'){ // '#' links a keyword with another keyword
                row_ind--;
                if (dualposition_flag)
                    copy_flag = 1; 
                else
                    copy_flag = 0;
                continue;
            }

            // cout << "i: "  << i << endl;
            // find a leaf node v and let the node i to be a child of node v
            int v = n;
            string childs_key(1, *q);
            while (nodes[v].childs.find(childs_key) != nodes[v].childs.end()){
                v = nodes[v].childs[childs_key];
                childs_key.append((++q), 1);
            }

            nodes[v].childs[childs_key] = i;

            if (copy_flag != 1){
                flag_count++;
                // FHE encrypt
                // nodes[i].bitmap[row] = 1;
                // bitmap_plaintext[row_ind] = 1;
                plain[i / ea_pointer->size()][i % ea_pointer->size()] = long2Poly(row_ind + 1);

                // cout << "before encrypt\n";
                // cout << "cipher.size() = " << cipher.size() << endl;
                // cout << "v / ea_pointer->size() = " << v / ea_pointer->size() << endl;
                // cout << "plain: " << plain << endl;
                // ea_pointer->encrypt(cipher[v / ea_pointer->size()], publicKey, plain);
                // ea_pointer->encrypt(cipher[0], publicKey, plain);
                // cout << "after encrypt\n";

                // bitmap_plaintext[row_ind] = 0;
                // plain[v % ea_pointer->size()] = 0;
                // int j = 0;
                // vector<NTL::ZZX>::iterator p = bitmap_plaintext.begin(); 
                // for (int j = 0; j < ptxt_num; j++) {
                //     vector<NTL::ZZX> temp_bitmap_paintext(p, p + ea_pointer->size());
                //     p += ea_pointer->size();
                //     ea_pointer->encrypt(*nodes[i].bitmap_cipher[j], publicKey, temp_bitmap_paintext);
                //     // ea_pointer->encrypt(*nodes[0].bitmap_cipher[j], publicKey, temp_bitmap_paintext);
                // }

                // encrypt plaintext to cipher 
                // std::vector<NTL::ZZX> plaintext;
                // plaintext.resize(ea_pointer->size(), NTL::ZZX(0));

                // int temp_ind = row_ind;
                // for (int y = 0; y < sizeof(int); y++){
                //     for (int z = 0; z < 8; z++){
                //         if (y*8 + z >= ea_pointer->size())
                //             break;
                //         plaintext[y*8 + z] = temp_ind;      
                //         temp_ind >>= 1;
                //     }
                // }
                
                // ea_pointer->encrypt(*nodes[i].cipher, publicKey, plaintext);


                // nodes[i].keyword.resize(records[row][attribute_index].size() + 2 * AES_BLOCK_SIZE);
                // RAND_bytes((unsigned char*)iv, AES_BLOCK_SIZE);
                // nodes[i].keyword.replace(0, AES_BLOCK_SIZE, (const char*)iv, AES_BLOCK_SIZE);
                // nodes[i].ciphertext_len = AES_encrypt((unsigned char *)&records[row][attribute_index][0], 
                //                                         records[row][attribute_index].size(), 
                //                                         aes_key, 
                //                                         iv, 
                //                                         (unsigned char *)&nodes[i].keyword[AES_BLOCK_SIZE]);
            }
        }
        // cout << cipher_num << endl;
        for (int i = 0; i < cipher_num; i++) {
            // cout << i << endl;
            // cout << "plain[i]:" << plain[i] << endl;
            ea_pointer->encrypt(cipher[i], publicKey, plain[i]);
        }
        // cout << "flag_count = " << flag_count << endl;
    }

    PositionHeap (const PositionHeap& p) {
        cout << "copy constructor\n";
    }

    PositionHeap& operator=(const PositionHeap& p) {
        cout << "operator = constructor\n";
        return *this;
    }

    ~PositionHeap() {
        // delete secretKey_pointer;
        // delete ea_pointer;
        // delete context_pointer;
        free(T);
    }

    void appendSubtree(FHESecKey *secretKey_pointer, 
                       int pos, 
                       vector<Ctxt>& matching_position, 
                       EncryptedArray *ea_pointer,
                    //    unsigned char *aes_key, 
                       int &time_count) {
        for (auto itr = nodes[pos].childs.begin(); itr != nodes[pos].childs.end(); itr++){
            // if (nodes[itr->second].keyword.size()){
            time_count++;
                // unsigned char output[nodes[itr->second].keyword.size()];
                // int plaintext_len;
                // plaintext_len = AES_decrypt((unsigned char*)&nodes[itr->second].keyword[AES_BLOCK_SIZE], 
                //                              nodes[itr->second].ciphertext_len, 
                //                              aes_key, 
                //                              (unsigned char*)&nodes[itr->second].keyword[0], 
                //                              output);
                // string plaintext((const char*)output, plaintext_len);
                // ret.push_back(plaintext);
            // temp_ret.push_back(nodes[itr->second].cipher);
            matching_position.push_back(cipher[itr->second / ea_pointer->size()]);

            // vector<NTL::ZZX> test_plaintext;
            // test_plaintext.resize(ea_pointer->size(), NTL::ZZX(0));
            // ea_pointer->decrypt(cipher[itr->second / ea_pointer->size()], *secretKey_pointer, test_plaintext);
            // cout << "line 311 test_plaintext: " << test_plaintext << endl;
            // cout << itr->second % ea_pointer->size() << endl;

            replicate(*ea_pointer, matching_position[matching_position.size() - 1], (itr->second) % ea_pointer->size());
            appendSubtree(secretKey_pointer, itr->second, matching_position, ea_pointer, time_count);
            // }
        }
    }
    
    void search (const char S[], 
                // unsigned char *aes_key, 
                int &time_count, 
                FHESecKey *secretKey_pointer, 
                EncryptedArray *ea_pointer,
                // vector<Ctxt *>& ret) {
                vector<Ctxt>& matching_position) {
        
        // vector<vector<Ctxt *>> temp_ret;
        int m = (int)strlen(S), depth = 0, v = n;
        string childs_key;
        // vector<int> matching_position;
        
        while (*S) {
            childs_key.append(S++, 1);
            if (nodes[v].childs.find(childs_key) == nodes[v].childs.end()){
                v = -1;
                break;
            }
            v = nodes[v].childs[childs_key];
            depth++;

            // if (nodes[v].keyword.size()){
            time_count++;
                // unsigned char output[nodes[v].keyword.size()];
                // int plaintext_len = AES_decrypt((unsigned char*)&nodes[v].keyword[AES_BLOCK_SIZE], nodes[v].ciphertext_len, 
                //                                 aes_key, (unsigned char*)&nodes[v].keyword[0], output);
                // string plaintext((const char*)output, plaintext_len); 
            // temp_ret.push_back(nodes[v].bitmap_cipher);
            // ret.push_back(nodes[v].cipher);

            matching_position.push_back(cipher[v / ea_pointer->size()]);

            // vector<NTL::ZZX> test_plaintext;
            // test_plaintext.resize(ea_pointer->size(), NTL::ZZX(0));
            // ea_pointer->decrypt(cipher[v / ea_pointer->size()], *secretKey_pointer, test_plaintext);
            // cout << "line 349 test_plaintext: " << test_plaintext << endl;
            // cout << v % ea_pointer->size() << endl;

            replicate(*ea_pointer, matching_position[matching_position.size() - 1], v % ea_pointer->size());
            // }
        }
        
        if (v != -1)
            appendSubtree(secretKey_pointer, v, matching_position, ea_pointer, time_count);


        // for (int i = 0; i < temp_ret.size(); i++){
        //     Ctxt *result = new Ctxt(*temp_ret[i]);
        //     ret.push_back(*result);
        // }

        // const FHEPubKey &publicKey = *secretKey_pointer;
        /* 
        if (temp_ret.size() == 0) 
            return;
        else{
            // for (int j = 0; j < ret.size(); j++)
            //     ret[j] = *temp_ret[0][j];
            ret[0] = *temp_ret[0];
        }

        // Todo: need to modfiy to multithread operations, which can reduce needed capcity
        // result = result | *ret[1] | *ret[2] | ...
        for (int i = 1; i < temp_ret.size(); i++){
            for (int j = 0; j < ret.size(); j++) {
                // note: Ctxt class hasn't operation* and operation+
                // ret[j] = ret[j] + (*temp_ret[j]) + ret[j] * (*temp_ret[j]);
                Ctxt temp = ret[j];
                ret[j] *= (*temp_ret[i][j]);
                ret[j] += (*temp_ret[i][j]);
                ret[j] += temp;
            }
        }

        cout << "ret[0] in search: " << ret[0] << endl;
        */

        // vector<NTL::ZZX> plaintextResult;

        // ea_pointer->decrypt(result, *secretKey_pointer, plaintextResult);

        // cout << "plaintextResult:" << plaintextResult << endl;

    }
};

#endif
