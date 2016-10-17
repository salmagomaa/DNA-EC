/**
 * A function that loads reads from file into memory
 * @param   String      orgFNm          The original reads file name
 * @return  int         0 for success & 1 for failure
 */
int getReads(string orgFNm) {

    ifstream fastq(orgFNm.c_str());

    if (fastq.is_open()) {

        int lnId = 1;
        string line;

        while (getline(fastq, line)) {
            if (lnId == 2) { //len
                if (readsCnt == 0) {
                    readsPtr = (struct Read *) malloc((readsCnt + 1) * sizeof (struct Read));
                } else {
                    readsPtr = (struct Read *) realloc(readsPtr, (readsCnt + 1) * sizeof (struct Read));
                }
                if (readsPtr == 0) {
                    printf("ERROR in allocating/reallocating readsPtr: Out of memory\n");
                    return 1;
                }
                readsPtr[readsCnt].len = line.length();
                readsPtr[readsCnt].seq = (char *) malloc((readsPtr[readsCnt].len + 1) * sizeof (char));
                if (readsPtr[readsCnt].seq == 0) {
                    printf("ERROR in allocating a seq: Out of memory\n");
                    return 1;
                }
                strcpy(readsPtr[readsCnt].seq, line.c_str());
            }else if(lnId == 4){
                readsPtr[readsCnt].qv = (char *) malloc((readsPtr[readsCnt].len + 1) * sizeof (char));
                if (readsPtr[readsCnt].qv == 0) {
                    printf("ERROR in allocating a qv: Out of memory\n");
                    return 1;
                }
                strcpy(readsPtr[readsCnt].qv, line.c_str());
            }
            if (lnId == seqLn) {
                readsCnt += 1;
                lnId = 0;
            }
            lnId++;
        }
        fastq.close();
    } else {
        cout << "Can't open the data file " << orgFNm << endl;
    }
    return 0;
}

/**
 * A function that set reads initial data (crrIdx)
 */
void setReadsData() {
    for (int i = 0; i < readsCnt; i++) {
        readsPtr[i].crr = (char *) malloc((readsPtr[i].len) * sizeof (char));
        readsPtr[i].crrFreq = (int *) malloc((readsPtr[i].len) * sizeof (int));
        if (readsPtr[i].crr == 0) {
            printf("ERROR in allocating a readsPtr[i].crr/readsPtr[i].crrFreq: Out of memory\n");
            return;
        }
        for (short int y = 0; y < readsPtr[i].len; y++) {
            readsPtr[i].crr[y] = CRR_NO;
            readsPtr[i].crrFreq[y] = 0;
        }
    }
}


/**
 * A function that dump reads' correction into a given file
 */
void dumpReads(string orgFNm, string crrFNm) {
    ifstream fastq(orgFNm.c_str());
    if (fastq.is_open()) {
        ofstream outFile(crrFNm.c_str(), ofstream::out);
        if (outFile.is_open()) {
            string line;
            int lnId = 1;
            int rIdx = -1;
            while (getline(fastq, line)) {
                if (lnId == 2) {//sequence
                    rIdx += 1;
                    line = readsPtr[rIdx].seq;
                }
                if (lnId == seqLn) {
                    lnId = 0;
                }
                lnId++;
                outFile << line << endl;
            }
            fastq.close();
            outFile.close();
            cout << "Reads are corrected and saved in " << crrFNm << " file!" << endl;
        } else {
            cout << "Can't create the corrections file " << crrFNm << endl;
        }
    } else {
        cout << "Can't re-open the data file " << orgFNm << endl;
    }
}

/**
 * A function that frees the allocated memory for reads
 */
void freeReads() {
    for (int i = 0; i < readsCnt; i++) {
        free(readsPtr[i].seq);
        free(readsPtr[i].crr);
        free(readsPtr[i].qv);
    }
    free(readsPtr);
    readsCnt = 0;
}