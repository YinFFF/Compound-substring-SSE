#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <map>            // std::map
#include <vector>         // std::vector 
#include <sys/time.h>     // gettimeofday
#include <stack>          // std::stack
#include <stdexcept>

#include "PHIndex.h"

#include <chrono>
#include <cmath>
#include <bitset>
#include <ctime>

#include <omp.h>

#include "NTL/ZZX.h"

#include "helib/FHE.h"
#include "helib/EncryptedArray.h"
#include "helib/binaryArith.h"
#include "helib/intraSlot.h"


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
    // if (B.capacity() < A.capacity()){
    //     Ctxt temp = B;
    //     B *= A;
    //     B += A; // B*A + A
    //     B += temp;  
    // } else {
    //     Ctxt temp = A;
    //     temp *= B;
    //     temp += B;
    //     temp += A;
    //     B = temp;
    // }
    Ctxt temp = B;
    B.multiplyBy(A);
    B += A;
    B += temp;
}

// B = B * A
void operationAND(Ctxt& A, Ctxt& B) {
    // if (B.capacity() < A.capacity())
    //     B *= A;
    // else {
    //     Ctxt temp = A;
    //     temp *= B;
    //     B = temp;
    // }
    B.multiplyBy(A);
}

// B[i] = A[i] + B[i] + c[i - 1], where c[i - 1] is carry-bit
// c[i] = (A[i] + c[i - 1]) * (B[i] + c[i - 1]) + c[i - 1]
// void operationAdd(vector<Ctxt>& A, vector<Ctxt>& B, const FHEPubKey &publicKey) {
//     Ctxt c(publicKey);
//     for (int i = 0; i < A.size(); i++) {
//         if (i > 0)
//             B[i] += c;
//         B[i] += A[i];

//         Ctxt temp1 = A[i];
//         temp1 += c;
//         Ctxt temp2 = B[i];
//         temp2 += c;
//         temp1 *= temp2;
//         c += temp1;
//     }
// }

// EqualTest for integers
void operationEQ(vector<Ctxt>& A, 
                vector<Ctxt>& B, 
                Ctxt& oneCiphertext, 
                Ctxt& result) {

    vector<Ctxt> temp_list;
    for (int i = 0; i < A.size() || i < B.size(); i++) {
        if (i < A.size() && i < B.size()) {
            Ctxt temp = A[i];
            temp += B[i];
            temp += oneCiphertext;
            temp_list.push_back(temp);
        } else if (i < A.size()) {
            Ctxt temp = A[i];
            temp += oneCiphertext;
            temp_list.push_back(temp);
        } else {
            Ctxt temp = B[i];
            temp += oneCiphertext;
            temp_list.push_back(temp);
        }
    }

    int right_ind = temp_list.size() - 1;
    while(right_ind) {
        int operation_count = (right_ind + 1) / 2;
        for (int i = 0; i < operation_count; ++i)
            // operationAND(temp_list[right_ind - i], temp_list[i]);
            temp_list[i].multiplyBy(temp_list[right_ind - i]);
        right_ind -= operation_count;
    }

    result = temp_list[0];

    // Ctxt oneCiphertext(publicKey);
    // std::vector<NTL::ZZX> onePlaintext(numSlots, NTL::ZZX(1));
    // ea.encrypt(oneCiphertext, publicKey, onePlaintext);
    // for (int i = 0; i < A.size() || i < B.size(); i++){
    //     if (i == 0) {
    //         result = oneCiphertext;
    //         result += A[i];
    //         result += B[i];
    //     } else if (i < A.size() && i < B.size()) {
    //         Ctxt temp = oneCiphertext;
    //         temp += A[i];
    //         temp += B[i];
    //         result *= temp;
    //     } else if (i < A.size()) {
    //         Ctxt temp = oneCiphertext;
    //         temp += A[i];
    //         result *= temp;
    //     } else {
    //         Ctxt temp = oneCiphertext;
    //         temp += B[i];
    //         result *= temp;
    //     }
    // }
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
    // long m = 255;
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

// \overline(t'_j) = EQ(\overline(id), \overline(id'))
// \overline(t'_j) = \overline(t_j) | \overline(t'_j)
void recursive_OR(const FHEPubKey &publicKey, 
                  EncryptedArray * ea_pointer, 
                  vector<vector<Ctxt>>& matching_id_sets, 
                  vector<vector<vector<Ctxt>>>& matching_pos_sets, 
                //   Ctxt& cur_cipher,
                  int cur_row_ind,
                  int cur_col_ind,
                  vector<Ctxt>& keyword_len_cipher,
                  Ctxt& oneCiphertext,
                  int row_ind,
                  int left_col_ind, 
                  int right_col_ind,
                  Ctxt& col_equal_flag){
    struct timeval time1, time2;
    float evaluate_time = 0;
    Ctxt temp_flag1(publicKey);
    Ctxt temp_flag2(publicKey);

    if (left_col_ind == right_col_ind){
        equalTest(temp_flag1, 
                matching_id_sets[cur_row_ind][cur_col_ind], 
                matching_id_sets[row_ind][left_col_ind], 
                ea_pointer->getDegree());

        vector<Ctxt> temp_flag2;
        CtPtrs_vectorCt pos_wrapper(temp_flag2);
        // cout << "   before addTwoNumbers\n";
        // gettimeofday(&time1,NULL);
        addTwoNumbers(
            pos_wrapper, 
            CtPtrs_vectorCt(keyword_len_cipher), 
            CtPtrs_vectorCt(matching_pos_sets[row_ind][left_col_ind]),
            false);
        // gettimeofday(&time2,NULL);
        // evaluate_time =  1000 * ((time2.tv_sec-time1.tv_sec)+((double)(time2.tv_usec-time1.tv_usec))/1000000);
        // printf("addTwoNumbers time: %f (ms)\n",  evaluate_time);

        // gettimeofday(&time1,NULL);
        operationEQ(temp_flag2, matching_pos_sets[cur_row_ind][cur_col_ind], oneCiphertext, col_equal_flag);
        // gettimeofday(&time2,NULL);
        // evaluate_time =  1000 * ((time2.tv_sec-time1.tv_sec)+((double)(time2.tv_usec-time1.tv_usec))/1000000);
        // printf("operationEQ time: %f (ms)\n",  evaluate_time);

        operationAND(temp_flag1, col_equal_flag);
        return;
    }

    int mid_col_ind = left_col_ind + (right_col_ind - left_col_ind) / 2;
    recursive_OR(
        publicKey, 
        ea_pointer, 
        matching_id_sets, 
        matching_pos_sets,
        cur_row_ind, 
        cur_col_ind, 
        keyword_len_cipher,
        oneCiphertext,
        row_ind, 
        left_col_ind, 
        mid_col_ind, 
        temp_flag1);
    recursive_OR(
        publicKey, 
        ea_pointer, 
        matching_id_sets, 
        matching_pos_sets,
        cur_row_ind, 
        cur_col_ind, 
        keyword_len_cipher,
        oneCiphertext,
        row_ind, 
        mid_col_ind + 1, 
        right_col_ind, 
        temp_flag2);

    operationOR(temp_flag1, temp_flag2);
    col_equal_flag = temp_flag2;
}

// filter for conjunctive logic
void filter_conjunction(const FHEPubKey &publicKey, 
            EncryptedArray * ea_pointer, 
            vector<vector<Ctxt>>& matching_id_sets, 
            vector<vector<vector<Ctxt>>>& matching_pos_sets, 
            vector<Ctxt> keyword_len_cipher,
            vector<Ctxt>& return_keywords){
            // int &small_size, 
            // int &largest_size){
    // find the small size of elment in matching_keywords

    // int small_ind = 0;
    // int largest_ind = 0;

    // make sure the smallest set and the largest set
    // for (int i = 0; i < matching_id_sets.size(); i++){
    //     if (matching_id_sets[i].size() < matching_id_sets[small_ind].size())
    //         small_ind = i;

    //     if (matching_id_sets[i].size() > matching_id_sets[largest_ind].size())
    //         largest_ind = i;
    // }
    
    // cout << "small_size = " <<  matching_id_sets[small_ind].size() << endl;
    // cout << "largest_size = " <<  matching_id_sets[largest_ind].size() << endl;
    // small_size = matching_id_sets[small_ind].size();
    // largest_size = matching_id_sets[largest_ind].size();


    std::vector<NTL::ZZX> onePlaintext(ea_pointer->size(), NTL::ZZX(1));
    Ctxt oneCiphertext(publicKey);
    ea_pointer->encrypt(oneCiphertext, publicKey, onePlaintext);

    // Ctxt temp_equal_flag(publicKey); 
    Ctxt t(publicKey); 
    Ctxt col_equal_flag(publicKey); 
    int equalTest_cnt = 0;
    // for (int j = 0; j < matching_id_sets[small_ind].size(); j++) {
    for (int j = 0; j < matching_id_sets[0].size(); j++) {
        cout  << "=========== filter the " << j << "th result in matching_id_sets[0]==============" << endl;
        t = oneCiphertext;
        // cout << "row_equal_flag.capacity():" << row_equal_flag.capacity() << endl;
        for (int row_ind = 0; row_ind < matching_id_sets.size(); row_ind++) {
            // col_equal_flag = zeroCiphertext;
            if (row_ind == 0 || matching_id_sets[row_ind].size() == 0)
                continue;
            recursive_OR(publicKey,
                        ea_pointer,
                        matching_id_sets,
                        matching_pos_sets,
                        // matching_id_sets[small_ind][j], 
                        0,
                        j,
                        keyword_len_cipher,
                        oneCiphertext,
                        row_ind,
                        0,
                        matching_id_sets[row_ind].size() - 1,
                        col_equal_flag);

            operationAND(col_equal_flag, t);
        }
        // Ctxt *new_keywords = new Ctxt(*matching_keywords[small_ind][j]);
        return_keywords.push_back(matching_id_sets[0][j]);

        // vector<NTL::ZZX> plaintextResult;
        // FHESecKey *secretKey_pointer = (FHESecKey *)&publicKey;
        // ea_pointer->decrypt(t, *secretKey_pointer, plaintextResult);
        // cout << "t:" << plaintextResult[0] << endl;
        // // return_keywords[return_keywords.size() - 1] *= row_equal_flag;
        // cout << "left capacity of t: " << t.capacity() << endl;
    }
}
    

// if DUALPOSITION = 1, scheme supports s1*s2 wildcard pattern query
#define DUALPOSITION 0
int TestNewSolution(char *file_name, int records_limit, int AttributesSize){

    // time calculation
    struct timeval time1, time2;
    float evaluate_time = 0;
    float evaluate_time1 = 0;
    float evaluate_time2 = 0;
    
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
    int pos_digit = 4;
    gettimeofday(&time1,NULL);
    vector<PositionHeap *> keywords_index;
    // for (int i = 0; i < AttributesSize; i++) {
    for (int i = 0; i < 1; i++) {
        PositionHeap *heap = new PositionHeap(string_set[i].c_str(),
                                              records, 
                                              i + 1,
                                              pos_digit,
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
    vector<vector<Ctxt>> matching_id_sets;
    vector<vector<vector<Ctxt>>> matching_pos_sets;
    vector<vector<float>> size_time;
    int small_size;
    int largest_size;

    // cout << "start query\n";
    for (int i = 0; i < 30; i++) {
        // cout << "matching_keywords.size(): " << matching_keywords.size() << endl;
        // cout << "========= i : " << i << " =========\n";
        int keyword_len = records[i][1].size();
        if (keyword_len <= 6)
            continue;
        vector<string> query_keywords;
        query_keywords.push_back(records[i][1].substr(0, (keyword_len - 1) / 2));
        query_keywords.push_back(records[i][1].substr((keyword_len - 1) / 2 + 1, (keyword_len - 1) / 2));
        cout << "query_keyword:" << query_keywords[0] << ", " << query_keywords[1] << endl;

        gettimeofday(&time1,NULL);
        for (int j = 0; j < query_keywords.size(); j++){            
            vector<Ctxt> temp_matching_id;
            vector<vector<Ctxt>> temp_matching_pos;
            keywords_index[0]->search(query_keywords[j].c_str(), 
            // keywords_index[0]->search(records[i][0 + 1].substr(0, ), 
                                        // aes_key, 
                                        time_count,
                                        secretKey_pointer,
                                        ea_pointer, 
                                        temp_matching_id, 
                                        temp_matching_pos);
            cout << "matching.size() = " << temp_matching_id.size() << endl;
            matching_id_sets.push_back(temp_matching_id);
            matching_pos_sets.push_back(temp_matching_pos);

            // print the matching records and pos
            // for (int i = 0; i < temp_matching_id.size(); i++) {

            //     vector<NTL::ZZX> plaintextResult(ea_pointer->size(), NTL::ZZX(0));

            //     ea_pointer->decrypt(temp_matching_id[i], *secretKey_pointer, plaintextResult);

            //     cout << "id:" << poly2Long(plaintextResult[0]) << endl;

            //     // cout << "temp_matching_pos[i].size():" << temp_matching_pos[i].size() << endl;
            //     vector<long> pos_result;
            //     CtPtrs_vectorCt matching_pos_wrapper(temp_matching_pos[i]);
            //     decryptBinaryNums(pos_result, matching_pos_wrapper, *secretKey_pointer, *ea_pointer);
            //     cout << "pos: " << pos_result[0] << endl;

            // }
        }
        gettimeofday(&time2, NULL);
        //msec
        evaluate_time1 = 1000*((time2.tv_sec-time1.tv_sec) + ((double)(time2.tv_usec-time1.tv_usec))/1000000);
        cout << "search time: " << evaluate_time1 <<  "(ms)" << endl;

        // encrypt keyword_len to keyword_len_cipher
        Ctxt scratch(publicKey);
        vector<Ctxt> keyword_len_cipher(pos_digit, scratch);
        for (int i = 0; i < pos_digit; i++) {
            vector<long> keyword_len_plain(ea_pointer->size());
            for (auto& slot : keyword_len_plain)
                slot = ((query_keywords[0].size() + 1) >> i) & 1;
            ea_pointer->encrypt(keyword_len_cipher[i], publicKey, keyword_len_plain);
        }
        gettimeofday(&time1, NULL);

        vector<Ctxt> result_set;
        filter_conjunction(publicKey, 
                ea_pointer, 
                matching_id_sets, 
                matching_pos_sets,
                keyword_len_cipher,
                result_set); 
                // small_size, 
                // largest_size);
        gettimeofday(&time2, NULL);

        evaluate_time2 = 1000*((time2.tv_sec-time1.tv_sec) + ((double)(time2.tv_usec-time1.tv_usec))/1000000);
        cout << "filter time: " << evaluate_time2 <<  "(ms)" << endl;
        
        vector<float> temp = {(float)matching_id_sets[0].size(), (float)matching_id_sets[1].size(), evaluate_time1, evaluate_time2};
        size_time.push_back(temp);

        //msec
        // evaluate_time = 1000*((time2.tv_sec-time1.tv_sec) + ((double)(time2.tv_usec-time1.tv_usec))/1000000);
        // cout << "filter time: " << evaluate_time << endl;

        matching_id_sets.clear();
        matching_pos_sets.clear();
        result_set.clear();
    }

    for (int i = 0; i < size_time.size(); i++)
        cout << size_time[i][0] << "\t" << size_time[i][1] << "\t"<< size_time[i][2] << "\t" << size_time[i][3] << endl;

    for (int i = 0; i < keywords_index.size(); i++)
        delete keywords_index[i];

    freeFHE(secretKey_pointer, context_pointer, ea_pointer);

    return 0;  
}


void IndexBuildForPreSolu(long numSlots, 
                long numEmpty, 
                long wordLength, 
                FHESecKey *secretKey_pointer,
                EncryptedArray* ea_pointer, 
                Ctxt &ciphertextAttr){
    const FHEPubKey &publicKey = *secretKey_pointer;

    char blankChar = 2;
    char wcChar = 3;
    char excludeChar = 4;

    // build the index
    // * wildcard is encoded as ASCII code 2
    // # wildcard is encoded as ASCII code 3
    std::string attrString = "spares"; // use the same attributes
    std::cout << std::endl << "Attribute String: " << attrString << std::endl;

    std::vector<NTL::ZZX> plaintextAttr; 
#if new_simd
    for(long i = 0; i < attrString.length(); i++)
        plaintextAttr.resize((i+1)*numWords, char2Poly(attrString[i]));
    for(long i = attrString.length(); i < wordLength; i++)
        plaintextAttr.resize((i+1)*numWords, plaintextAttr[i] = char2Poly(blankChar));
#else
    plaintextAttr.resize(numSlots-numEmpty, NTL::ZZX(0));
    for (unsigned long i = 0; i < attrString.length(); i++) 
        plaintextAttr[i] = char2Poly(attrString[i]);
    for (unsigned long i = attrString.length(); i < wordLength; i++)
        plaintextAttr[i] = char2Poly(blankChar);
    for (unsigned long i = wordLength; i < numSlots - numEmpty; i++)
        plaintextAttr[i] = plaintextAttr[i % wordLength];    
#endif
    plaintextAttr.resize(numSlots, NTL::ZZX(0));

    // Ctxt ciphertextAttr(publicKey);
    ea_pointer->encrypt(ciphertextAttr, publicKey, plaintextAttr);

    return;
}


int TestPreSolution(char *file_name, int records_limit, int AttributesSize){

    // Timers
    auto startTime = std::chrono::high_resolution_clock::now();
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> timeTaken = endTime - startTime;
    float totalTime = 0;

    // initial FHE
    FHESecKey *secretKey_pointer = NULL;
    FHEcontext *context_pointer = NULL;
    EncryptedArray *ea_pointer = NULL;
    InitializeFHE(secretKey_pointer, context_pointer, ea_pointer);

    long numSlots = ea_pointer->size(); // the number of points in a ciphertext or plaintext
    long plaintextDegree = ea_pointer->getDegree();
    const FHEPubKey &publicKey = *secretKey_pointer;

    char blankChar = 2;
    char wcChar = 3;
    char excludeChar = 4;

    long wordLength = 10; 
    long queryLength;      
     
    long numWords = floor(numSlots / wordLength); 
    long numEmpty = numSlots % wordLength;
    long conjSize = 3;//atoi(argv[4]);

    Ctxt ciphertextAttr(publicKey);
    IndexBuildForPreSolu(numSlots, numEmpty, wordLength, secretKey_pointer, ea_pointer, ciphertextAttr); 
    
    // "$" is the wildcard character symbol
    std::string queryString = "sp";
    queryString += wcChar;
    queryString += excludeChar;
    queryString += "c";
    queryString += "e";
    queryLength = queryString.length();
    std::cout << "Query Pattern: " << queryString << std::endl;

    std::vector<NTL::ZZX> plaintextQuery, plaintextE;
    // std::vector<NTL::ZZX> plaintextConjunction((unsigned long) numSlots, long2Poly(rand() % (1L << 7)));

    plaintextQuery.resize(numSlots-numEmpty, NTL::ZZX(0));
    plaintextE.resize(numSlots-numEmpty, NTL::ZZX(0));
    int counter = 0;
    for (unsigned long i = 0; i < queryString.length(); i++) {
        if (queryString[i] == wcChar) {
            plaintextQuery[counter] = char2Poly(excludeChar);
            plaintextE[counter] = 1;
            counter++;
        } else if (queryString[i] == excludeChar) {
            plaintextQuery[counter] = char2Poly(queryString[i + 1]);
            plaintextE[counter] = 1;
            i += 1;
            counter++;
        } else {
            plaintextQuery[counter] = char2Poly(queryString[i]);
            counter++;
        }
    }
    for (unsigned long i = counter; i < wordLength; i++) {
        plaintextQuery[i] = char2Poly(wcChar);
        plaintextE[i] = 1;
    }
    for (unsigned long i = wordLength; i < numSlots - numEmpty; i++) {
        plaintextQuery[i] = plaintextQuery[i % wordLength];
        plaintextE[i] = plaintextE[i % wordLength];
    }

    plaintextQuery.resize(numSlots, NTL::ZZX(0));
    plaintextE.resize(numSlots, NTL::ZZX(0));

    std::vector<NTL::ZZX> plaintextResult((unsigned long) numSlots - numEmpty, NTL::ZZX(1));
    plaintextResult.resize(numSlots, NTL::ZZX(0));

    std::clog << "wordLength:" << wordLength << ", " << "queryLength: " << queryLength << ", " << "numWords: " << numWords << ", ";

    // Initialize and encrypt query ciphertexts
    Ctxt ciphertextQuery(publicKey);
    Ctxt tempCiphertext(publicKey);
    Ctxt ciphertextE(publicKey);
    Ctxt ciphertextResult(publicKey);
    Ctxt conjResult(publicKey);
//        Ctxt conjQuery(publicKey);
//        Ctxt conjCtxt(publicKey);

    ea_pointer->encrypt(ciphertextQuery, publicKey, plaintextQuery);
    ea_pointer->encrypt(ciphertextE, publicKey, plaintextE);
//        ea.encrypt(conjCtxt, publicKey, plaintextConjunction);
//        ea.encrypt(conjQuery, publicKey, plaintextConjunction);

    Ctxt oneCiphertext(publicKey);
    Ctxt zeroCiphertext(publicKey);

    std::vector<NTL::ZZX> onePlaintext(numSlots, NTL::ZZX(1));

    ea_pointer->encrypt(oneCiphertext, publicKey, onePlaintext);
    zeroCiphertext = oneCiphertext;
    zeroCiphertext -= oneCiphertext;

    // ciphertextWs is a set of \overline{A}
    // 文章的 algorithm 1 中有一个次数为 wordlength 的循环， 这里通过 vector 的形式统一运算
    std::vector<Ctxt> ciphertextWs, ciphertextRs, ciphertextSs, ciphertextDs;
    for (unsigned long i = 0; i < wordLength; i++) {
        ciphertextWs.push_back(ciphertextAttr);
        ciphertextRs.push_back(ciphertextAttr);
        ciphertextSs.push_back(ciphertextAttr);
    }
    
    // ciphertextDs is a set of \overline{D} 
    for (unsigned long i = 0; i < wordLength; i++) {
        if (i <= wordLength - queryLength) {
            ciphertextDs.push_back(oneCiphertext);
        } else {
            ciphertextDs.push_back(zeroCiphertext);
        }
    }

    // For compound conjunction queries
    std::vector<Ctxt> ciphertextConj;
    std::vector<Ctxt> ciphertextConjResult;
    for (unsigned long i = 0; i < conjSize; i++) {
        ciphertextConj.push_back(ciphertextAttr);
        ciphertextConjResult.push_back(ciphertextAttr);
    }

    std::cout << "Encryption Done!" << std::endl;

    NTL::ZZX selectPoly;
    NTL::ZZX queryMask, finalMask;
    makeMask(selectPoly, 1, 1, wordLength, *ea_pointer);
    makeMask(finalMask, wordLength, 1, wordLength, *ea_pointer);

    // Step 1: Shift the attributes
    //     Remnants of experiment to pack differently, all characters of the same slot first
    //     But the shift time is much longer than packing word by word

    //     startTime = std::chrono::high_resolution_clock::now();
    // #pragma omp parallel for
    //     for(unsigned long i = 1; i < ciphertextWs.size(); i++) {
    //             ea.shift(ciphertextWs[i],-i*numWords);
    //     }
    //     endTime = std::chrono::high_resolution_clock::now();
    //     timeTaken = endTime-startTime;
    //     totalTime += timeTaken.count();
    //     std::cout << "Pre-compute Time: " << timeTaken.count() << std::endl;

    startTime = std::chrono::high_resolution_clock::now();
    for (unsigned long i = 1; i < ciphertextWs.size(); i++) {
        simdShift(ciphertextWs[i], ciphertextAttr, -i, wordLength, *ea_pointer);
    }
    endTime = std::chrono::high_resolution_clock::now();
    timeTaken = endTime - startTime;
    // totalTime += timeTaken.count();
    std::cout << "Pre-compute Time: " << timeTaken.count() << std::endl;

    // Step 2: Test if the characters are the same
    startTime = std::chrono::high_resolution_clock::now();
    // X = A - W
    for (unsigned long i = 0; i < ciphertextWs.size() + conjSize; i++) {
        if (i < ciphertextWs.size()) {
            ciphertextWs[i] += ciphertextQuery;
        } else {
            ciphertextConj[i - ciphertextWs.size()] += ciphertextAttr;
        }
    }
    endTime = std::chrono::high_resolution_clock::now();
    timeTaken = endTime - startTime;
    totalTime += timeTaken.count();
    std::cout << "XOR Time: " << timeTaken.count() << std::endl;

    startTime = std::chrono::high_resolution_clock::now();
    // Y = EQ(x, 0)
    for (unsigned long i = 0; i < ciphertextWs.size() + conjSize; i++) {
        if (i < ciphertextWs.size()) {
            // ciphertextRs is a set of \overline{Y}
            equalTest(ciphertextRs[i], zeroCiphertext, ciphertextWs[i], plaintextDegree);
        } else {
            equalTest(ciphertextConjResult[i - ciphertextWs.size()], zeroCiphertext,
                      ciphertextConj[i - ciphertextWs.size()], plaintextDegree);
        }
        // S <- Y_i + E
        // for(unsigned long i = 0; i < ciphertextRs.size(); i++) {
        // ciphertextRs is \overline{S}
        if (i < ciphertextWs.size()) {
            ciphertextRs[i] += ciphertextE;
        }
    }
    endTime = std::chrono::high_resolution_clock::now();
    timeTaken = endTime - startTime;
    totalTime += timeTaken.count();
    std::cout << "Level left: " << ciphertextRs[1].capacity() << ", Equality Check + eMask Time: " << timeTaken.count() << std::endl;
    std::clog << timeTaken.count() << ", ";

    // Step 3: Combine results of character tests per shift
    startTime = std::chrono::high_resolution_clock::now();
    // z_i = \prod_{s_i,j}
    for (unsigned long i = 0; i < ciphertextRs.size() + conjSize; i++) {
        if (i < ciphertextRs.size()) {
            oneTotalProduct(ciphertextSs[i], ciphertextRs[i], wordLength, *ea_pointer);
        } else if (conjSize > 0) {
            oneTotalProduct(ciphertextConj[i - ciphertextRs.size()], ciphertextConjResult[i - ciphertextRs.size()],
                            wordLength, *ea_pointer);
        }
    }
    endTime = std::chrono::high_resolution_clock::now();
    timeTaken = endTime - startTime;
    totalTime += timeTaken.count();
    std::cout << "Level left: " << ciphertextSs[1].capacity() << ", Product of Equalities: " << timeTaken.count() << std::endl;
    std::clog << timeTaken.count() << ", ";

    // Step 4: Combine results from testing every possible shift
    startTime = std::chrono::high_resolution_clock::now();
    // 1 + z_i * d_i
    // ciphertextSs is a set of \overline{z_i}
    for (unsigned long i = 0; i < ciphertextSs.size(); i++) {
        ciphertextSs[i].multiplyBy(ciphertextDs[i]);
        ciphertextSs[i].addConstant(selectPoly);
    }
    endTime = std::chrono::high_resolution_clock::now();
    timeTaken = endTime - startTime;
    totalTime += timeTaken.count();
    std::cout << "Level left: " << ciphertextSs[1].capacity() << ", Disjunction Prep Time: " << timeTaken.count() << std::endl;
    std::clog << timeTaken.count() << ", ";

    startTime = std::chrono::high_resolution_clock::now();
    // ciphertextResult = 1 + \sigma(1 + z * d)
    treeMult(ciphertextResult, ciphertextSs);
    ciphertextResult.addConstant(selectPoly);

    // \sigma * \sigma * \belta_{i,j} + 1
    if (conjSize > 0) {
        treeMult(conjResult, ciphertextConj);
        ciphertextResult.multiplyBy(conjResult);
    }
    // ciphertextResult.multByConstant(finalMask);
    endTime = std::chrono::high_resolution_clock::now();
    timeTaken = endTime - startTime;
    totalTime += timeTaken.count();
    std::cout << "Level left: " << ciphertextResult.capacity() << ", Disjunction Time: " << timeTaken.count() << std::endl;
    std::clog << timeTaken.count() << ", ";

    std::clog << totalTime << ", ";

    startTime = std::chrono::high_resolution_clock::now();
    tempCiphertext = ciphertextResult;
    
    for (unsigned long i = 1; i < wordLength; i++) {
        ciphertextWs[i] = ciphertextResult;
        ea_pointer->shift(ciphertextWs[i], i);
    }
    for (unsigned long i = 1; i < wordLength; i++) {
        tempCiphertext += ciphertextWs[i];
    }
    tempCiphertext.multiplyBy(ciphertextAttr);
    // 此时， 如果 ciphertextResult 中第一个位置是1，那么 tempCiphertext 是 “sparce”
    endTime = std::chrono::high_resolution_clock::now();
    timeTaken = endTime - startTime;
    // totalTime += timeTaken.count();
    std::cout << "Level left: " << tempCiphertext.capacity() << ", Fill + Selection Time: " << timeTaken.count() << std::endl;
    // std::clog << timeTaken.count() << ", ";

    std::cout << "Total Time: " << totalTime << std::endl;

    std::cout << "ciphertextResult: " << std::endl;
    ea_pointer->decrypt(ciphertextResult, *secretKey_pointer, plaintextResult);
    for (unsigned long i = 0; i < wordLength; i++) {
        std::cout << poly2Long(plaintextResult[i]) << ", ";
    }
    std::cout << std::endl;

    std::cout << "tempciphertext: " << std::endl;
    ea_pointer->decrypt(tempCiphertext, *secretKey_pointer, plaintextResult);
    for (unsigned long i = 0; i < wordLength; i++) {
        std::cout << poly2Char(plaintextResult[i]) << ", ";
    }
    std::cout << std::endl;

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
            TestNewSolution(file_name, record_num2, 2);
            // cout << "=============== TestPreSolution =================\n";
            // TestInvertedIndex(file_name, 5000, AttrSize);
            // cout << "=============== TestSuffixTreeSolution =================\n";
            // TestSuffixTreeSolution(file_name, record_num, AttrSize);
            // getchar();
        }
    }
    return 0;
}


