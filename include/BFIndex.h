#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H

#include <stdio.h>
#include <iostream>
#include "bloom_filter.hpp"

using namespace std;

class test_class
{

    int a;
    char *c;
    int b;
    public:
        test_class(int a, int b){
            this->a = a;
            this->b = b;
            c = new char[10];
        }
        ~test_class(){
            delete[] c;
        }
};



class BFIndex { 
    bloom_parameters parameters;
    bloom_filter *filter;
    string ciphertext_keyword;
    int ciphertext_len;
    
public:
    BFIndex(string &plaintext_keyword, unsigned char *aes_key){
        
        struct timeval time1,time2;

        // How many elements roughly do we expect to insert?
        parameters.projected_element_count = 1000;
 
        // Maximum tolerable false positive probability? (0,1)
        parameters.false_positive_probability = 0.0001; // 1 in 10000
 
        // Simple randomizer (optional)
        parameters.random_seed = 0xA5A5A5A5;
 
        if (!parameters){
            std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
        }

        parameters.compute_optimal_parameters();

        filter = new bloom_filter(parameters); 
       
        for (int j = 0; j < plaintext_keyword.size(); j++)
            for (int z = j; z < plaintext_keyword.size(); z++)
                filter->insert(plaintext_keyword.substr(j, z - j + 1));

        ciphertext_keyword = plaintext_keyword;
        
        ciphertext_keyword.resize(plaintext_keyword.size() + AES_BLOCK_SIZE);
        
        unsigned char iv[AES_BLOCK_SIZE];
        
        RAND_bytes((unsigned char*)iv, AES_BLOCK_SIZE);
        
        ciphertext_keyword.replace(0, AES_BLOCK_SIZE, (const char*)iv, AES_BLOCK_SIZE);

        ciphertext_len = AES_encrypt((unsigned char *)&plaintext_keyword[0], 
                                                    plaintext_keyword.size(), 
                                                    aes_key, 
                                                    iv, 
                                                    (unsigned char *)&ciphertext_keyword[AES_BLOCK_SIZE]);
                                                    
    }
    
    ~BFIndex(){
        delete filter;
    }

    string search(string &search_keyword, unsigned char *aes_key) {
        string result;
        
        if (filter -> contains(search_keyword)){
            unsigned char output[ciphertext_keyword.size()];
            
            int plaintext_len = AES_decrypt((unsigned char*)&ciphertext_keyword[AES_BLOCK_SIZE], ciphertext_len, 
                        aes_key, (unsigned char*)&ciphertext_keyword[0], output);
            
            result.assign((const char*)output, plaintext_len);   
        }
        
        return result;
    }        
};

#endif
