#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <map>            // std::map
#include <vector>         // std::vector 
#include <sys/time.h>     // gettimeofday
#include <stack>          // std::stack
//#include <openssl/hmac.h>
//#include <openssl/aes.h>
//#include <openssl/rand.h>
#include <stdexcept>

#include "PHIndex.h"
// #include "BFIndex.h"
// #include "AES.h"
// #include "BMIndex.h"
// #include "SuffixTree.h"

#include <chrono>
#include <cmath>
#include <bitset>
#include <ctime>

#include <omp.h>

#include "NTL/ZZX.h"

#include "helib/FHE.h"
#include "helib/EncryptedArray.h"

using namespace std;

// plaintext 长度必须为 numSlots
#define new_simd 0

NTL::ZZX long2Poly(long query) {
    NTL::ZZX queryPoly;
    queryPoly.SetLength(64);
    std::bitset<64> queryBits = query;

    for (int i = 0; i < queryBits.size(); i++) {
        NTL::SetCoeff(queryPoly, i, queryBits[i]);
    }
    queryPoly.normalize();

    return queryPoly;
}

NTL::ZZX char2Poly(char character) {
    int charCode = character;
    NTL::ZZX resultPoly;

    resultPoly = long2Poly((long) charCode);

    return resultPoly;
}

long poly2Long(NTL::ZZX result) {
    long resultLong = 0;
    for (int i = 0; i <= deg(result); i++) {
        resultLong += (1L << i) * (NTL::coeff(result, i) == NTL::ZZ(1));
    }
    return resultLong;
}

char poly2Char(NTL::ZZX result) {
    char character;

    character = (int) poly2Long(result);

    return character;
}

void fastPower(Ctxt &dataCtxt, long degree) {
    // Taken from eqtesting.cpp so that there are fewer includes
    if (degree == 1) return;
    Ctxt orig = dataCtxt;
    long k = NTL::NumBits(degree);
    long e = 1;
    for (long i = k - 2; i >= 0; i--) {
        Ctxt tmp1 = dataCtxt;
        tmp1.smartAutomorph(1L << e); // 1L << e computes 2^e
        dataCtxt.multiplyBy(tmp1);
        e = 2 * e;
        if (NTL::bit(degree, i)) {
            dataCtxt.smartAutomorph(2);
            dataCtxt.multiplyBy(orig);
            e += 1;
        }
    }
}

void equalTest(Ctxt &resultCtxt, const Ctxt &queryCtxt, const Ctxt &dataCtxt,
               long degree) {
    Ctxt tempCtxt = dataCtxt;
    // tempCtxt = dataCtxt;
    tempCtxt -= queryCtxt;
    fastPower(tempCtxt, degree);
    tempCtxt.negate();
    tempCtxt.addConstant(NTL::ZZ(1));
    resultCtxt = tempCtxt;
}

void treeMultHelper(Ctxt &resultCtxt, std::vector<Ctxt> &currentLayer, std::vector<Ctxt> &nextLayer) {
    unsigned long previousSize = currentLayer.size();
    if (previousSize == 0) {
        return;
    } else if (previousSize == 1) {
        resultCtxt = currentLayer[0];
        return;
    }
    nextLayer.resize((previousSize / 2 + previousSize % 2), resultCtxt);
#pragma omp parallel for
    for (unsigned long i = 0; i < previousSize / 2; i++) {
        currentLayer[2 * i].multiplyBy(currentLayer[2 * i + 1]);
        nextLayer[i] = currentLayer[2 * i];
    }

    if (previousSize % 2 == 1) {
        nextLayer[nextLayer.size() - 1] = (currentLayer[previousSize - 1]);
    }
    currentLayer.clear();
    treeMultHelper(resultCtxt, nextLayer, currentLayer);
}

void treeMult(Ctxt &resultCtxt, const std::vector<Ctxt> &ctxtVec) {
    if (ctxtVec.size() > 1) {
        std::vector<Ctxt> currentLayer, nextLayer;
        currentLayer = ctxtVec;
        nextLayer.clear();
        treeMultHelper(resultCtxt, currentLayer, nextLayer);
    } else {
//        std::cout << "Only 1 Ciphertext; No Multiplication Done." << std::endl;
        resultCtxt = ctxtVec[0];
    }
}

void oneTotalProduct(Ctxt &resultCtxt, const Ctxt &dataCtxt, const long wordLength, const EncryptedArray &ea) {

    long numWords = floor(ea.size() / wordLength);
    resultCtxt = dataCtxt;
    if (wordLength == 1) {
        return;
    }
    long shiftAmt = 1;

    // auto startTime = std::chrono::high_resolution_clock::now();
    // auto endTime = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> timeTaken = endTime - startTime;

    while (shiftAmt < wordLength) {
        // startTime = std::chrono::high_resolution_clock::now();
        Ctxt tempCtxt = resultCtxt;
#if new_simd
        ea.shift(tempCtxt, (-shiftAmt*numWords));
#else
         ea.shift(tempCtxt, (-shiftAmt));
#endif
        resultCtxt.multiplyBy(tempCtxt); // ctxt = ctxt * (ctxt << "shiftAmt")

        shiftAmt = 2 * shiftAmt;
        // endTime = std::chrono::high_resolution_clock::now();
        // timeTaken = endTime-startTime;
        // std::cout << shiftAmt << ", Time Taken: " << timeTaken.count() << std::endl;
    }
}

void makeMask(NTL::ZZX &maskPoly, const long shiftAmt, const bool invertSelection, const long wordLength,
              const EncryptedArray &ea) {
    std::vector<NTL::ZZX> maskVec, oneMaskVec;
    if (invertSelection) {
        maskVec.assign(shiftAmt, NTL::ZZX(1));
        oneMaskVec.assign(wordLength - shiftAmt, NTL::ZZX(0));
    } else {
        maskVec.assign(shiftAmt, NTL::ZZX(0));
        oneMaskVec.assign(wordLength - shiftAmt, NTL::ZZX(1));
    }
    maskVec.insert(maskVec.end(), oneMaskVec.begin(), oneMaskVec.end());

    std::vector<NTL::ZZX> fullMaskVec = maskVec;
    for (unsigned long i = 2 * wordLength; i < ea.size(); i += wordLength) {
        fullMaskVec.insert(fullMaskVec.end(), maskVec.begin(), maskVec.end());
    }

    fullMaskVec.resize(ea.size(), NTL::ZZX(0));
    ea.encode(maskPoly, fullMaskVec);
}

void simdShift(Ctxt &ciphertextResult, Ctxt &ciphertextData, const long shiftAmt, const long wordLength,
               const EncryptedArray &ea) {
    Ctxt tempCiphertext = ciphertextData;
    if (shiftAmt > 0) {
        NTL::ZZX maskPoly;
        makeMask(maskPoly, shiftAmt, 0, wordLength, ea);

        ea.shift(tempCiphertext, shiftAmt);
        tempCiphertext.multByConstant(maskPoly);
    } else if (shiftAmt < 0) {
        NTL::ZZX maskPoly;
        makeMask(maskPoly, -shiftAmt, 0, wordLength, ea);

        tempCiphertext.multByConstant(maskPoly);
        ea.shift(tempCiphertext, shiftAmt);
    }
    ciphertextResult = tempCiphertext;
}

// B = B * A + A + B
void operationOR(Ctxt& A, Ctxt& B) {
    Ctxt temp = B;
    B *= A;
    B += A; // B*A + A
    B += temp;  
}

// B = B * A
void operationAND(Ctxt& A, Ctxt& B) {
    B *= A;
}

/*
	srcStr:原始string
	Keyword:分割符
	destVect:分割后放入容器中
 */ 
int splitStringToVect(const string & srcStr, vector<string> & destVect, const string & Keyword)  
{  
	string Tempstr;			//to store current substring
	int Pos = 0;			
	int Templen;		
	while(1)
	{
		Templen = srcStr.find(Keyword, Pos) - Pos;		
		if (Templen < 0)								
			break;
        else if (Templen > 0){
            Tempstr = srcStr.substr(Pos, Templen);			
			destVect.push_back(Tempstr);
        }
		Pos += Templen + Keyword.length();				
	}
    if(Pos < srcStr.length())
        destVect.push_back(srcStr.substr(Pos, srcStr.length() - Pos)); //The last substring
    return destVect.size();  
}


void readKeywords(vector<vector<string> > &attributes_list, char* file_name, int records_size, int attributes_size){
    string line;
    ifstream readfile(file_name);
    int count = 0;
    if(readfile.is_open()){
        while(count != records_size && getline(readfile, line, '\n')){
            if(line.find('#') != -1 or line.find('_') != -1) throw invalid_argument("find a '#' or '_' in file\n");
            if (line[line.size() - 1] == '\r')
                line.erase(line.end()  - 1, line.end());
            vector<string> keywords;
            splitStringToVect(line, keywords, " ");
            if (keywords.size() < attributes_size + 1)
                continue;
            attributes_list.push_back(keywords);
            count ++;
        }
        readfile.close();
    }
    else
        printf("file open error\n");
    return;
}

int keywords_to_str(vector<vector<string>> &records, int attributes_size, vector<string> &string_set, bool flag){

    for (int j = 1; j < attributes_size + 1; j++){
        int total_len = 0;
        string keywords_string;
        
        for(int i = 0; i < records.size(); i++) {
            if (flag)
                total_len += 2 * (records[i][j].size());
            else
                total_len += records[i][j].size();
        }

        if (flag)
            keywords_string.resize(total_len + 2 * records.size());
        else
            keywords_string.resize(total_len + records.size());
        
        int current_position = 0;
        for(int i = 0; i < records.size(); i++){
            keywords_string.replace(current_position, records[i][j].size(), records[i][j]);
            current_position += records[i][j].size();
            // copy attribute
            if (flag){
                keywords_string.replace(current_position, 1, "_");
                current_position ++;
                keywords_string.replace(current_position, records[i][j].size(), records[i][j]);
                current_position += records[i][j].size();
                keywords_string.replace(current_position, 1, "#");
                current_position ++;
            } else{
                keywords_string.replace(current_position, 1, "#");
                current_position ++;
            }
        }
        
        string_set.push_back(keywords_string);
    }
    return 1;
}

struct time_count{
    float total_time;
    int test_num;
};


int InitializeFHE(FHESecKey*& secretKey_pointer, FHEcontext*& context_pointer, EncryptedArray*& ea_pointer){

    // FHE initialize
    long p = 2;
    long r = 1;
    long m = 18631;//21607;//20485; //18631;//32767; //16383;//32767; // 8191;//32767;
    int cap = 500;
    // long L = 31; //31
    // long m = FindM(80, 200, 2, p, 0, 0, 0);

    context_pointer = new FHEcontext(m, p, r);
    FHEcontext& context = *context_pointer;
    // cout << "context addr: " << &context << endl;
    buildModChain(context, cap);

    NTL::ZZX F = context.alMod.getFactorsOverZZ()[0];

    secretKey_pointer = new FHESecKey(context);

    FHESecKey &secretKey = *secretKey_pointer;
    const FHEPubKey &publicKey = secretKey;
    secretKey.GenSecKey(64);
    addFrbMatrices(secretKey);
    addBSGS1DMatrices(secretKey);
    ea_pointer = new EncryptedArray(context, F);
    // cout << "ea_pointer->getcontext addr: " << &ea_pointer->getContext() << endl;

    long numSlots = ea_pointer->size();
    long plaintextDegree = ea_pointer->getDegree();

    long wordLength = 10;//atoi(argv[3]);
    long queryLength = 17;
    long numWords = floor(numSlots / wordLength);
    long numEmpty = numSlots % wordLength;

    IndexSet allPrimes(0, context.numPrimes() - 1);
    std::clog << "m=" << m << ", " << "cap=" << cap << ", " 
    << context_pointer->logOfProduct(context_pointer->ctxtPrimes) / log(2.0) 
    << ", " <<  context_pointer->logOfProduct(allPrimes) / log(2.0) 
    << ", " << "p=" << p << ", " << plaintextDegree
    << ", " << "numSlots=" << numSlots 
    << ", " << "security= " << context.securityLevel() << endl;
    return 0;
}

void freeFHE(FHESecKey* secretKey_pointer, FHEcontext* context_pointer, EncryptedArray* ea_pointer) {
    delete secretKey_pointer;
    delete context_pointer;
    delete ea_pointer;
}

void inter_to_ptxt(EncryptedArray* ea_pointer, int inter, vector<NTL::ZZX> &ptxt) {
    // ptxt.resize(ea_pointer->size(), NTL::ZZX(0));
    int temp_ind = inter;
    for (int i = 0; i < sizeof(int); i++){
        for (int j = 0; j < 8; j++){
            if (i*8 + j >= ea_pointer->size())
                break;
            ptxt[i*8 + j] = temp_ind;      
            temp_ind >>= 1;
        }
    }
}

void ptxt_to_inter(EncryptedArray* ea_pointer, vector<NTL::ZZX> &ptxt, int &inter) {
    long temp;
    inter = 0;
    for(int i = sizeof(int) - 1; i >= 0; i--){
        for (int j = 7; j >= 0; j--){
            if (i*8 + j >= ea_pointer->size())
                continue;
            inter <<= 1;
            temp = poly2Long(ptxt[i*8 + j]);
            inter |= temp;
        }
    }
}

void recursive_OR(const FHEPubKey &publicKey, 
                  EncryptedArray * ea_pointer, 
                  vector<vector<Ctxt>>& matching_positions, 
                  Ctxt& cur_cipher,
                  int row_ind,
                  int left_col_ind, 
                  int right_col_ind,
                  Ctxt& col_equal_flag){
    // struct timeval time1, time2;
    // float evaluate_time = 0;
    if (left_col_ind == right_col_ind){
        // gettimeofday(&time1,NULL);
        equalTest(col_equal_flag, cur_cipher, matching_positions[row_ind][left_col_ind], ea_pointer->getDegree());
        // gettimeofday(&time2,NULL);
        // evaluate_time =  1000 * ((time2.tv_sec-time1.tv_sec)+((double)(time2.tv_usec-time1.tv_usec))/1000000);
        // printf("equalTest: %f (ms)\n",  evaluate_time);
        return;
    }

    Ctxt temp_col_equal_flag(publicKey);
    int mid_col_ind = left_col_ind + (right_col_ind - left_col_ind) / 2;
    recursive_OR(publicKey, ea_pointer, matching_positions, cur_cipher, row_ind, left_col_ind, mid_col_ind, temp_col_equal_flag);
    recursive_OR(publicKey, ea_pointer, matching_positions, cur_cipher, row_ind, mid_col_ind + 1, right_col_ind, col_equal_flag);

    // gettimeofday(&time1,NULL);
    operationOR(temp_col_equal_flag, col_equal_flag);
    // gettimeofday(&time2,NULL);
    // evaluate_time =  1000 * ((time2.tv_sec-time1.tv_sec)+((double)(time2.tv_usec-time1.tv_usec))/1000000);
    // printf("operationOR: %f (ms)\n",  evaluate_time);
}

void filter(const FHEPubKey &publicKey, 
            EncryptedArray * ea_pointer, 
            vector<vector<Ctxt>>& matching_positions, 
            vector<Ctxt>& return_keywords,
            int &small_size, 
            int &largest_size){
    // find the small size of elment in matching_keywords

    int small_ind = 0;
    int largest_ind = 0;

    // make sure the smallest set and the largest set
    for (int i = 0; i < matching_positions.size(); i++){
        if (matching_positions[i].size() < matching_positions[small_ind].size())
            small_ind = i;

        if (matching_positions[i].size() > matching_positions[largest_ind].size())
            largest_ind = i;
    }
    
    cout << "small_size = " <<  matching_positions[small_ind].size() << endl;
    cout << "largest_size = " <<  matching_positions[largest_ind].size() << endl;
    small_size = matching_positions[small_ind].size();
    largest_size = matching_positions[largest_ind].size();


    std::vector<NTL::ZZX> onePlaintext(ea_pointer->size(), NTL::ZZX(1));
    Ctxt oneCiphertext(publicKey);
    ea_pointer->encrypt(oneCiphertext, publicKey, onePlaintext);

    // Ctxt temp_equal_flag(publicKey); 
    Ctxt row_equal_flag(publicKey); 
    Ctxt col_equal_flag(publicKey); 

    int equalTest_cnt = 0;
    int AND_cnt = 0;
    for (int j = 0; j < matching_positions[small_ind].size(); j++) {
        cout  << j << "th result" << endl;
        row_equal_flag = oneCiphertext;
        // cout << "row_equal_flag.capacity():" << row_equal_flag.capacity() << endl;
        for (int row_ind = 0; row_ind < matching_positions.size(); row_ind++) {
            // col_equal_flag = zeroCiphertext;
            if (row_ind == small_ind || matching_positions[row_ind].size() == 0)
                continue;

            recursive_OR(publicKey, ea_pointer, matching_positions, matching_positions[small_ind][j], row_ind, 0, matching_positions[row_ind].size() - 1, col_equal_flag);
            // cout << "after recursive_OR:"  << col_equal_flag.capacity() << endl;

            operationAND(col_equal_flag, row_equal_flag);
            // cout << "operationAND " << AND_cnt << ": " << col_equal_flag.capacity() << "\t" << row_equal_flag.capacity() << endl;
            AND_cnt++;
        }
        // Ctxt *new_keywords = new Ctxt(*matching_keywords[small_ind][j]);
        return_keywords.push_back(matching_positions[small_ind][j]);
        // return_keywords[return_keywords.size() - 1] *= row_equal_flag;
        cout << "left capacity: " << row_equal_flag.capacity() << endl;
    }
}
    

// if DUALPOSITION = 1, scheme supports s1*s2 wildcard pattern query
#define DUALPOSITION 0
int TestNewSolution(char *file_name, int records_limit, int AttributesSize){

    // time calculation
    struct timeval time1, time2;
    float evaluate_time = 0;
    
    // read keywords to keywords_list from the file
    vector<vector<string>> records;
    readKeywords(records, file_name, records_limit, AttributesSize); 

    // transfer records to a set of string
    vector<string> string_set;
    keywords_to_str(records, AttributesSize, string_set, DUALPOSITION);

    // initial FHE
    FHESecKey *secretKey_pointer = NULL;
    FHEcontext *context_pointer = NULL;
    EncryptedArray *ea_pointer = NULL;
    InitializeFHE(secretKey_pointer, context_pointer, ea_pointer);
    const FHEPubKey &publicKey = *secretKey_pointer;

    if (records.size() >= pow(2, ea_pointer->getDegree())){
        cout << "Degree is too small\n";
        return 1;
    }
        
    // Build a set of position heap
    gettimeofday(&time1,NULL);
    vector<PositionHeap *> keywords_index;
    for (int i = 0; i < AttributesSize; i++){
        PositionHeap *heap = new PositionHeap(string_set[i].c_str(),
                                            //   aes_key, 
                                              records, 
                                              i + 1,
                                              secretKey_pointer, 
                                              ea_pointer, 
                                              context_pointer,
                                              DUALPOSITION);
        keywords_index.push_back(heap);
    }
    gettimeofday(&time2,NULL);
   
    //msec
    evaluate_time =  1000 * ((time2.tv_sec-time1.tv_sec)+((double)(time2.tv_usec-time1.tv_usec))/1000000);
    printf("Positioinheap build time: %f (ms)\n",  evaluate_time);

    // vector<string>().swap(string_set);
    // for (int i = 0; i < records.size(); i++){
    //     vector<string>().swap(records[i]);
    // }

    // Input query keywords
    // vector<string> query_keywords;
    // for (int j = 0; j < AttributesSize; j++){
    //     string query_keyword;
    //     cout << "input keyword " << j << ":";
    //     getline(cin, query_keyword);
    //     query_keywords.push_back(query_keyword);
    // }

    int ptxt_num = ceil((float) records.size() / ea_pointer->size());
    int time_count = 0;
    vector<vector<Ctxt>> one_attribute_result;

    vector<vector<float>> size_time;
    int small_size;
    int largest_size;

    cout << "start query\n";
    for(int i = 100; i < 120; i++) {
        // cout << "matching_keywords.size(): " << matching_keywords.size() << endl;
        cout << "========= i : " << i << " =========\n";
        gettimeofday(&time1,NULL);
        for (int j = 0; j < AttributesSize; j++){            
            vector<Ctxt> temp_matching_position;
            // keywords_index[j]->search(query_keywords[j].c_str(), 
            keywords_index[j]->search(records[i][j + 1].c_str(), 
                                        // aes_key, 
                                        time_count,
                                        secretKey_pointer,
                                        ea_pointer, 
                                        temp_matching_position);
            // cout << "temp_matching_position.size() = " << temp_matching_position.size() << endl;
            one_attribute_result.push_back(temp_matching_position);
        }

        vector<Ctxt> multi_attributes_result;
        filter(publicKey, ea_pointer, one_attribute_result, multi_attributes_result, small_size, largest_size);
        gettimeofday(&time2, NULL);

        //msec
        evaluate_time = 1000*((time2.tv_sec-time1.tv_sec) + ((double)(time2.tv_usec-time1.tv_usec))/1000000);
        cout << "search time: " << evaluate_time << endl;
        
        vector<float> temp = {(float)small_size, (float)largest_size, evaluate_time};
        size_time.push_back(temp);

        //msec
        // evaluate_time = 1000*((time2.tv_sec-time1.tv_sec) + ((double)(time2.tv_usec-time1.tv_usec))/1000000);
        // cout << "filter time: " << evaluate_time << endl;

        one_attribute_result.clear();
        multi_attributes_result.clear();
    }

    for (int i = 0; i < size_time.size(); i++)
        cout << size_time[i][0] << "\t" << size_time[i][1] << "\t"<< size_time[i][2] << endl;

    for (int i = 0; i < keywords_index.size(); i++)
        delete keywords_index[i];
    freeFHE(secretKey_pointer, context_pointer, ea_pointer);

    return 0;  
}


int main(int argc, char * argv[])
{
    if (argc != 3){
	cout << "parameter error\n";
        return 1;
    }
    // the format must be unix style (ending style: '\n') 
    // char file_name[30] = "./Testfile/test";
    char file_name[50] = "./Testfile/records_with_5_attributes";
    // int record_num = 2000;
    int AttrSize2 = atoi(argv[1]);
    int record_num2 = atoi(argv[2]);
    for (int AttrSize = 2; AttrSize <= 2; AttrSize += 1){
        cout << "========AttributeSize = " << AttrSize2 << "==========\n";
        for (int record_num = 500; record_num <= 500; record_num += 500) {
            cout << "=============== TestNewSolution =================\n"; 
            cout << "records_num = " << record_num2 << endl;
            TestNewSolution(file_name, record_num2, AttrSize2);
            // cout << "=============== TestPreSolution =================\n";
            // TestInvertedIndex(file_name, 5000, AttrSize);
            // cout << "=============== TestSuffixTreeSolution =================\n";
            // TestSuffixTreeSolution(file_name, record_num, AttrSize);
            // getchar();
        }
    }
    return 0;
}


