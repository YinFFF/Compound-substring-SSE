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
    int ciphertext_len;

    Node(int bitmap_size,  const FHEPubKey& publicKey){
        ciphertext_len = 0;
    }

    Node (const Node& p){
        childs = p.childs;
        keyword = p.keyword;
        // bitmap = p.bitmap;
        ciphertext_len = p.ciphertext_len;
    }

    Node (Node&& p) {
        // cout << "3" << endl;
        childs = p.childs;
        keyword = p.keyword;
        // bitmap = p.bitmap;
        ciphertext_len = p.ciphertext_len;
    }

    Node& operator=(const Node& p){
        // cout << "4" << endl;
        childs = p.childs;
        keyword = p.keyword;
        // bitmap = p.bitmap;
        ciphertext_len = p.ciphertext_len;
        return *this;
    }

    Node& operator=(Node&& p){
        // cout << "5" << endl;
        if (this != &p){
            childs = p.childs;
            keyword = p.keyword;
            // bitmap = p.bitmap;
            ciphertext_len = p.ciphertext_len;
        }
        return *this;
    }

    ~Node(){
    }

};


class PositionHeap {
    char *T;
    int n;
    
public:
    vector<Node> nodes;
    vector<Ctxt> cipher; // cipher is used to store all the encrypted id stored in this position heap;
    vector<vector<Ctxt>> pos_cipher; // pos is used to store related positions of nodes

    PositionHeap(const char T[], 
                //  unsigned char *aes_key, 
                 vector<vector<string>> &records, 
                 int attribute_index,
                 int pos_digit,
                 FHESecKey *secretKey_pointer,
                 EncryptedArray *ea_pointer,
                 FHEcontext *context_pointer,
                 int dualposition_flag) {

        const FHEPubKey &publicKey = *secretKey_pointer;
        // unsigned char iv[AES_BLOCK_SIZE];
        int row_ind = records.size();

        // int id_division = row_ind;
        // int id_remainder = 0;
        // int id_digit = 0;

        // while(id_division) {
        //     id_remainder = id_division % 2;
        //     id_division = id_division / 2;
        //     id_digit++; 
        // }

        // int pos_digit = 8;

        // int ptxt_num = ceil((float)row_ind / ea_pointer->size());
        int copy_flag;
        n = (int)strlen(T);
        int cipher_num = ceil((float) (n + 1) / ea_pointer->size());
        this->T = strdup(T);
        cout << "n:" << n << endl;
        cout << "T: " << T << endl;


        nodes.resize(n+1, Node(ea_pointer->size(), publicKey));
        cipher.resize(cipher_num, Ctxt(publicKey));
        vector<Ctxt> temp_pos_cipher(pos_digit, Ctxt(publicKey));
        pos_cipher.resize(cipher_num, temp_pos_cipher);


        // vector<NTL::ZZX> bitmap_plaintext(ea_pointer->size() * ptxt_num, NTL::ZZX(0));
        // vector<vector<NTL::ZZX>> plain;
        // for (int i = 0; i < cipher_num; i++) {
        //     vector<NTL::ZZX> temp_plain;
        //     temp_plain.resize(ea_pointer->size(), NTL::ZZX(0));
        //     plain.push_back(temp_plain);
        // }

        vector<NTL::ZZX> id_plain;
        id_plain.resize(ea_pointer->size(), NTL::ZZX(0));
        // pos_plain represents the current 
        vector<vector<NTL::ZZX>> pos_plain; 
        for (int i = 0; i < pos_digit; i++) {
            vector<NTL::ZZX> temp_pos_plain;
            temp_pos_plain.resize(ea_pointer->size(), NTL::ZZX(0));
            pos_plain.push_back(temp_pos_plain);
        }

        int flag_count = 0;
        int pos = 0; // the distance between current position and the end character in a keyword
        
        for (int i = n; --i >= 0; pos++) {
            const char *q = &T[i];
            if ((*q) == '_'){  // '_' links a keyword with its cope
                copy_flag = 0; //  
                pos = 0;
                continue;
            }else if ((*q) == '#'){ // '#' links a keyword with another keyword
                row_ind--;
                pos = 0;
                if (dualposition_flag)
                    copy_flag = 1; 
                else
                    copy_flag = 0;
                continue;
            }

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
                // plain[i / ea_pointer->size()][i % ea_pointer->size()] = long2Poly(row_ind + 1);
                id_plain[i % ea_pointer->size()] = long2Poly(row_ind + 1);

                // put pos to pos_plain
                int pos_division = pos;
                for(int j = 0; j < pos_digit && pos_division; j++){
                    pos_plain[j][i % ea_pointer->size()] = pos_division % 2;
                    pos_division = pos_division / 2;
                }
            }

            if (i % ea_pointer->size() == 0){
                ea_pointer->encrypt(cipher[i / ea_pointer->size()], publicKey, id_plain);

                // encrypt pos_plain to pos_cipher[i / ea_pointer->size()]
                for(int j = 0; j < pos_digit; j++){
                    ea_pointer->encrypt(pos_cipher[i / ea_pointer->size()][j], publicKey, pos_plain[j]);
                    for(int z = 0; z < ea_pointer->size(); z++)
                        pos_plain[j][z] = 0;
                }
            }

        }

        // FHE encrypt
        // for (int i = 0; i < cipher_num; i++) {
        //     ea_pointer->encrypt(cipher[i], publicKey, plain[i]);
        // }
    }

    PositionHeap (const PositionHeap& p) {
        cout << "copy constructor\n";
    }

    PositionHeap& operator=(const PositionHeap& p) {
        cout << "operator = constructor\n";
        return *this;
    }

    ~PositionHeap() {
        free(T);
    }

    void appendSubtree(FHESecKey *secretKey_pointer, 
                       int pos, 
                       vector<Ctxt>& matching_id, 
                       vector<vector<Ctxt>>& matching_position, 
                       EncryptedArray *ea_pointer,
                    //    unsigned char *aes_key, 
                       int &time_count) {
        for (auto itr = nodes[pos].childs.begin(); itr != nodes[pos].childs.end(); itr++){
            // if (nodes[itr->second].keyword.size()){
            time_count++;
            matching_id.push_back(cipher[itr->second / ea_pointer->size()]);
            matching_position.push_back(pos_cipher[itr->second / ea_pointer->size()]);

            // replicate(*ea_pointer, matching_id[matching_id.size() - 1], (itr->second) % ea_pointer->size());
            ea_pointer->shift(matching_id[matching_id.size() - 1], - ((itr->second) % ea_pointer->size()));
            for (int i = 0; i < matching_position[matching_position.size() - 1].size(); i++)
                // replicate(*ea_pointer, matching_position[matching_position.size() - 1][i], (itr->second) % ea_pointer->size());
                ea_pointer->shift(matching_position[matching_position.size() - 1][i], - ((itr->second) % ea_pointer->size()));
            appendSubtree(secretKey_pointer, itr->second, matching_id, matching_position, ea_pointer, time_count);
            // }
        }
    }
    
    void search (const char S[], 
                // unsigned char *aes_key, 
                int &time_count, 
                FHESecKey *secretKey_pointer, 
                EncryptedArray *ea_pointer,
                // vector<Ctxt *>& ret) {
                vector<Ctxt>& matching_id,
                vector<vector<Ctxt>>& matching_position) {
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

            matching_id.push_back(cipher[v / ea_pointer->size()]);
            matching_position.push_back(pos_cipher[v / ea_pointer->size()]);

            // vector<NTL::ZZX> test_plaintext;
            // test_plaintext.resize(ea_pointer->size(), NTL::ZZX(0));
            // ea_pointer->decrypt(cipher[v / ea_pointer->size()], *secretKey_pointer, test_plaintext);
            // cout << "line 349 test_plaintext: " << test_plaintext << endl;
            // cout << v % ea_pointer->size() << endl;

            // replicate(*ea_pointer, matching_id[matching_id.size() - 1], v % ea_pointer->size());
            ea_pointer->shift(matching_id[matching_id.size() - 1], - (v % ea_pointer->size()));
            for (int i = 0; i < matching_position[matching_position.size() - 1].size(); i++)
                // replicate(*ea_pointer, matching_position[matching_position.size() - 1][i], v % ea_pointer->size());
                ea_pointer->shift(matching_position[matching_position.size() - 1][i], - (v % ea_pointer->size()));
        }
        
        if (v != -1)
            appendSubtree(secretKey_pointer, v, matching_id, matching_position, ea_pointer, time_count);
    }
};

#endif
