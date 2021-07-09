# Pre-requisites

The program has been successfully built and tested on CentOS 6.6. It uses OpenSSL library to implement the entailed cryptographic operations. The build commands in CentOS are as follows.

## OpenSSL

```sh
$ wget http://www.openssl.org/source/openssl-1.1.1.tar.gz
$ tar -xvf openssl-1.1.1.tar.gz
$ cd openssl-1.1.1/
$ ./configure
$ make
$ make install
```



## Getting the code
The code is available *via* git:

```sh
 $ git clone git@github.com:YinFFF/Substring-keyword-SSE.git
```

# Usage

This project includes two source files: AES.cpp and SubstringQueryScheme.cpp: the former contains interface functions for AES, and the latter contains the implement of the test functions.

## Build and run

```sh
$ cd Substring-keyword-SSE/
$ ./make
$ ./SubstringQueryScheme

```