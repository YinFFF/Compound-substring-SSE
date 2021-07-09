#ifndef BITMAP_H
#define BITMAP_H

#include <stdio.h>
#include <iostream>
#include <map>
#include "AES.h"
#include <string>

using namespace std;

class BMIndex { 
public:
    int id;
    int n; // the nubmer of bytes in a row of the bitmap
    map<string, unsigned char *> bitmap;

    BMIndex(vector<vector<string>> &records, int attribute_index, unsigned char* records_key){
        n = records.size() / 8 + 1;
        for (id = 0; id < records.size(); id++){
            unsigned char * Hmac_output1 = NULL;  
	        unsigned int len_of_Hmac_output1 = 0; 
//            HmacEncode("sha256", (const char*)records_key, 32, records[id][attribute_index].c_str(), 
//                records[id][attribute_index].size(), Hmac_output1, len_of_Hmac_output1);
//            string K((char *)Hmac_output1, len_of_Hmac_output1);
//            free(Hmac_output1);
            string K(records[id][attribute_index]);
            if (bitmap.find(K) == bitmap.end()){
                bitmap[K] = new unsigned char[n];
                memset(bitmap[K], 0, n * sizeof(char));
            }
            unsigned char *p = bitmap[K] + id / 8;
            *p = (*p) | (1 << (7 - id % 8));
        }
    }
    
    ~BMIndex(){
        for (auto itr = bitmap.begin(); itr != bitmap.end(); itr++)
            delete itr->second;
    }

    int search(string &search_keyword, unsigned char* records_key, unsigned char* output) {
        unsigned char * Hmac_output = NULL;  
        unsigned int len_of_Hmac_output = 0; 
//        HmacEncode("sha256", (const char*)records_key, 32, search_keyword.c_str(), search_keyword.size(), 
//		           Hmac_output, len_of_Hmac_output);
//        string K((char *)Hmac_output, len_of_Hmac_output);
        string K(search_keyword);
        if(bitmap.find(K) != bitmap.end()){
            memcpy(output, bitmap[K], n);
        }
        return 1;
    }        
};



class InvertedIndex{
    map<string, int> index;
public:
    InvertedIndex(vector<vector<string>> &records, int attribute_index, unsigned char* records_key){
        map<string, int> keyword_count;
        for (int i = 0; i < records.size(); i++){
            if (keyword_count.find(records[i][attribute_index]) == keyword_count.end())
                keyword_count[records[i][attribute_index]] = 1;
            else
                keyword_count[records[i][attribute_index]] += 1;        
            index[records[i][attribute_index] + to_string(static_cast<long 
            long>(keyword_count[records[i][attribute_index]]))] = i;
       }

//        cout << attribute_index << endl;
//        for (auto itr = index.begin(); itr != index.end(); itr++)
//            cout << itr->first << ":" << itr->second << endl;
        
    }

    vector<int> search(string &search_keyword, unsigned char* records_key){
        vector<int> result;
        int keyword_count = 1;
        while (index.find(search_keyword + to_string(static_cast<long long>(keyword_count))) != index.end()){
            result.push_back(index[search_keyword + to_string(static_cast<long long>(keyword_count))]);
            keyword_count++;
        }
        return result;
    }
};
#endif
