#ifndef SUFFIXTREE_H
#define SUFFIXTREE_H

#include <algorithm>
#include <exception>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <deque>
#include <map> 
#include <functional>
#include "AES.h"

using namespace std;


class SuffixTree {
    class Node {
    private:
        int ciphertext_len;
        string keyword;
    public:
        // size_t edge_hash;
        vector<int> ch;    // vector of child nodes
        int left_index;
        int right_index;

        Node() {
            keyword = "";
            left_index = -1;
            right_index = -1;
            ciphertext_len = 0;
        }

        Node(int left_index, int right_index, 
             std::initializer_list<int> children, 
             int pos, const string& suf, const string& keyword, 
             unsigned char *aes_key) {
            ch.insert(ch.end(), children);
            this->left_index = left_index;
            this->right_index = right_index;
            // string hash_string = suf.substr(0, left_index);
            // edge_hash = hash<string>{}(hash_string);

            if(keyword.size()){
                unsigned char iv[AES_BLOCK_SIZE];
                RAND_bytes((unsigned char*)iv, AES_BLOCK_SIZE);
                this->keyword.resize(keyword.size() + 2*AES_BLOCK_SIZE);
                this->keyword.replace(0, AES_BLOCK_SIZE, (const char*)iv, AES_BLOCK_SIZE);
                this->ciphertext_len = AES_encrypt((unsigned char *)&(keyword[0]), 
                                                        keyword.size(), 
                                                        aes_key, 
                                                        iv, 
                                                        (unsigned char *)&(this->keyword[AES_BLOCK_SIZE]));

            }else{
                this->ciphertext_len = 0;
            }
            
        }

        string getKeyword(unsigned char * aes_key){
            // decrypt keyword
            if (this->ciphertext_len){
                // unsigned char output[100];
                // cout << "keyword.size(): "<< this->keyword.size() << endl;
                // cout << "cipherte_len:" << this->ciphertext_len << endl;
                // cout << "enc_flag: " << enc_flag << endl;

                // cout << this->keyword.size() << " ";
                // printf("this->keyword: ");
                // for (int z = 0; z < this->keyword.size(); z++)
                //     printf("%0X", this->keyword[z]);
                // printf("\n");
                unsigned char output[keyword.size()];

                int plaintext_len = AES_decrypt((unsigned char*)&(this->keyword[AES_BLOCK_SIZE]), 
                                                this->ciphertext_len,
                                                aes_key, 
                                                (unsigned char*)&(this->keyword[0]), 
                                                output);
                string plaintext((const char*)output, plaintext_len); 
                return plaintext;
            }else
                return "";
        }
    };

public:
    vector<Node> nodes;
    
    SuffixTree(string & str, unsigned char *aes_key, vector<vector<string>> &records, int attribute_index) {
        int keyword_index = records.size();
        if (str[str.length() - 1] != '$')
            str.append("$");
        nodes.push_back(Node{});
        // subs.push_back("");
        for (int pos = str.length() - 2; pos >= 0; pos--) {
            if (str[pos] == '#') {
                keyword_index--;
                continue;
            }
            addSuffix(str, pos, aes_key, records[keyword_index][attribute_index]);
        }
    }

    void addSuffix(const std::string & suf, int pos, unsigned char *aes_key, const string & keyword) {
        int n = 0, x2, n2;

        for (int i = pos; i < suf.length(); ) {
            x2 = 0;
            while (true) {
                vector<int> children = nodes[n].ch;
                if (x2 == children.size()) {
                    // no matching child, remainder of suf becomes new node
                    n2 = nodes.size();
                    nodes.push_back(Node(i, suf.size() - 1, {}, pos, suf, keyword, aes_key));
                    nodes[n].ch.push_back(n2);
                    return;
                }
                n2 = children[x2];
                if (suf[nodes[n2].left_index] == suf[i]) {
                    break;
                }
                x2++;
            }
            // find prefix of remaining suffix in common with child
            int sub2_left_index = nodes[n2].left_index;
            int sub2_right_index = nodes[n2].right_index;
            size_t j = 0;
            while (j < sub2_right_index - sub2_left_index + 1) {
                if (suf[i + j] != suf[sub2_left_index + j]) {
                    // split n2
                    int n3 = n2;
                    // new node for the part in common
                    n2 = nodes.size();
                    nodes.push_back(Node(sub2_left_index, sub2_left_index + j - 1, {n3}, -1, suf, "", aes_key));
                    nodes[n3].left_index = sub2_left_index + j;
                    nodes[n3].right_index = suf.size() - 1;
                    nodes[n].ch[x2] = n2;
                    break; // continue down the tree
                }
                j++;
            }
            i += j; // advance past part in common
            n = n2; // continue down the tree
        }
    }

    void recursion_print(int n, const std::string & pre){
        auto children = nodes[n].ch;
        if (children.size() == 0) {
            // std::cout << "- " << nodes[n].left_index <<', ' << nodes[n].right_index << '\n';
            std::cout << "- " << nodes[n].left_index  << '\n';
            // std::cout << "- " << subs[n] << '\n';
            return;
        }

        std::cout << "+ " << nodes[n].left_index << '\n';
        // std::cout << "+ " << nodes[n].left_index <<', ' << nodes[n].right_index << '\n';
        // std::cout << "+ " << subs[n] << '\n';

        auto it = children.begin();
        if (it != children.end()) do {
            if ( (it + 1) == children.end()) break;
            std::cout << pre << "+-";
            recursion_print(*it, pre + "| ");
            it++;
        } while (true);

        std::cout << pre << "+-";
        recursion_print(children[children.size() - 1], pre + "  ");
    }

    void visualize() {
        if (nodes.size() == 0) {
            std::cout << "<empty>\n";
            return;
        }
        recursion_print(0, "");
    }

    void appendSubtree(int cur_node_index, unsigned char *aes_key, vector<string> &ret) {
        static int count = 0;
        if (nodes[cur_node_index].ch.size() == 0) {
            string keyword = nodes[cur_node_index].getKeyword(aes_key);
            if (keyword.size())
                ret.push_back(keyword);
            return;
        }
        
        for (auto p = nodes[cur_node_index].ch.begin(); p != nodes[cur_node_index].ch.end(); p++){
            appendSubtree(*p, aes_key, ret);
        }
    }    

    vector<string> search(const char S[], const string& suf, unsigned char *aes_key) {
        vector<string> ret;
        int m = strlen(S);
        int cur_node_index = 0, S_index = 0;

        // find the node, whose edges matching the queried string S
        while (S_index < m) {
            auto chs = nodes[cur_node_index].ch;
            int j;

            for (j = 0; j < chs.size(); j++) {
                // cout << "chs[j]:" << chs[j] << endl;
                if (suf[nodes[chs[j]].left_index] == S[S_index]) {
                    int sub_len = min(nodes[chs[j]].right_index - nodes[chs[j]].left_index + 1, m - S_index);
                    if (suf.compare(nodes[chs[j]].left_index, sub_len, &S[S_index], sub_len) == 0){
                        S_index += sub_len;
                        cur_node_index = chs[j];
                        break;
                    }
                    else
                        return ret;
                }
            }

            if (j == chs.size()){
                return ret;
            }
        }


        appendSubtree(cur_node_index, aes_key, ret);
        return ret;
    }


};
#endif
