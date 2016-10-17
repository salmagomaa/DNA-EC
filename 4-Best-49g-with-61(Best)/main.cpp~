#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include <stdint.h>
#include <tr1/unordered_map>

using namespace std;

ofstream meta("meta", ofstream::out);

#include "definitions.cpp"
#include "reads.cpp"
#include "kmers.cpp"
#include "correct.cpp"

int main() {
    struct timeval startedAt, finishedAt, crrStartedAt, crrFinishedAt;
    long mtime, seconds, useconds;
    gettimeofday(&startedAt, NULL);

    string orgFNm = "/home/ubuntu/SRR088759.fastq";
    short int levels =5;
    short int k = 25;
    checkDel = false;
    seqLn = 4;
    if (getReads(orgFNm) == 0) {

        meta << "Reads count: " << readsCnt << endl;

        gettimeofday(&crrStartedAt, NULL);

        setReadsData();
        meta << "Reads data are set successfully!" << endl;

        for (short int i = 0; i < levels; i++) {
            meta << "Level #" << i << endl;
            buildValidKmersHashing(k);
            meta << "Kmers count: " << kmersCnt << endl;

            correctMismatchedKmers(k);
            meta << "Mismatched kmers corrected successfully" << endl;
            
            correctSolidKmers(k);
            meta << "Solid kmers corrected successfully" << endl;

            kmersInvalidCnt = 0;
            for (int kk = 0; kk < kmersCnt; kk++) {
                if (kmersPtr[kk].validByKIdx == -1) {
                    kmersInvalidCnt += 1;
                }
            }
            meta << "kmersInvalidCnt = " << kmersInvalidCnt << endl;

            updateReadsByCorrections();
            meta << "Reads are updated with corrections successfully" << endl;
            
            if (i < (levels - 1)) {
                freeKmers();
                meta << "Kmers are free!" << endl;
            }
            meta << "*************" << endl;
        }
        
        gettimeofday(&crrFinishedAt, NULL);
        seconds = crrFinishedAt.tv_sec - crrStartedAt.tv_sec;
        useconds = crrFinishedAt.tv_usec - crrStartedAt.tv_usec;
        mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
        meta << "Correction consumed time in ms:" << mtime << endl;

        string crrFNm = "CorrectedReads";
        dumpReads(orgFNm, crrFNm);
        orgFNm = crrFNm;

        freeKmers();
        meta << "Kmers are free!" << endl;

        freeReads();
        meta << "Reads are free!" << endl;
    }
    gettimeofday(&finishedAt, NULL);
    seconds = finishedAt.tv_sec - startedAt.tv_sec;
    useconds = finishedAt.tv_usec - startedAt.tv_usec;
    mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
    meta << "Total consumed time in ms:" << mtime << endl;

    meta.close();

    return 0;
}
//////////////////////////////////////
