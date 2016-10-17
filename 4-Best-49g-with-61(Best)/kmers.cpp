/**
 *  Coral: short reads error correction with multiple alignments
 *  Copyright (C) 2011 Leena Salmela <leena.salmela@cs.helsinki.fi>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * A function that adds a read to a kmer's reads list
 * @param int   kIdx    The kmer idx within the kmers array
 * @param bool  isNew   A flag that represents whether given kmer is a new or already exists before
 * @param int   rdIdx   The read idx within the reads array
 * @param int   pos     The first idx of the kmer within the read's sequence
 * @param bool  isRev   A flag that represents whether the found kmer is a reverse complement or not 
 */
void addReadToKmer(int kIdx, bool isNew, int rdIdx, short int pos, bool isRev) {
    if (isNew) {
        kmersCnt += 1;
        if (kmersCnt == 1) {
            kmersPtr = (struct Kmer *) malloc(kmersCnt * sizeof (struct Kmer));
        } else {
            kmersPtr = (struct Kmer *) realloc(kmersPtr, kmersCnt * sizeof (struct Kmer));
        }
        if (kmersPtr == 0) {
            printf("ERROR in allocating/reallocating kmersPtr: Out of memory\n");
            return;
        }
        kmersPtr[kIdx].rCnt = 0;
        kmersPtr[kIdx].revRCnt = 0;
        kmersPtr[kIdx].validByKIdx = -1;
    }
    if (isRev) {
        kmersPtr[kIdx].revRCnt += 1;
        if (kmersPtr[kIdx].revRCnt == 1) {
            kmersPtr[kIdx].revRIds = (int *) malloc(sizeof (int));
            kmersPtr[kIdx].revStPos = (short int *) malloc(sizeof (short int));
        } else {
            kmersPtr[kIdx].revRIds = (int *) realloc(kmersPtr[kIdx].revRIds, kmersPtr[kIdx].revRCnt * sizeof (int));
            kmersPtr[kIdx].revStPos = (short int *) realloc(kmersPtr[kIdx].revStPos, kmersPtr[kIdx].revRCnt * sizeof (short int));
        }
        if (kmersPtr[kIdx].revRIds == 0 || kmersPtr[kIdx].revStPos == 0) {
            printf("ERROR in allocating/reallocating a revRIds/revStPos: Out of memory\n");
            return;
        }
        kmersPtr[kIdx].revRIds[kmersPtr[kIdx].revRCnt - 1] = rdIdx;
        kmersPtr[kIdx].revStPos[kmersPtr[kIdx].revRCnt - 1] = pos;
    } else {
        kmersPtr[kIdx].rCnt += 1;
        if (kmersPtr[kIdx].rCnt == 1) {
            kmersPtr[kIdx].rIds = (int *) malloc(sizeof (int));
            kmersPtr[kIdx].stPos = (short int *) malloc(sizeof (short int));
        } else {
            kmersPtr[kIdx].rIds = (int *) realloc(kmersPtr[kIdx].rIds, kmersPtr[kIdx].rCnt * sizeof (int));
            kmersPtr[kIdx].stPos = (short int *) realloc(kmersPtr[kIdx].stPos, kmersPtr[kIdx].rCnt * sizeof (short int));
        }
        if (kmersPtr[kIdx].rIds == 0 || kmersPtr[kIdx].stPos == 0) {
            printf("ERROR in allocating/reallocating a rIds/stPos: Out of memory\n");
            return;
        }
        kmersPtr[kIdx].rIds[kmersPtr[kIdx].rCnt - 1] = rdIdx;
        kmersPtr[kIdx].stPos[kmersPtr[kIdx].rCnt - 1] = pos;
    }
}

int hashKmer(char *str, int len, int rIdx, bool excludeN, short int k, bool getIdxOnly, bool *getIsRev) {

    int kIdx = -1;
    int bits = 2;
    int lastN = -1;
    int pos = -1;
    uint64_t addingGram = 0;
    uint64_t gram = 0;
    uint64_t gram_reverse = 0;
    uint64_t mask = ((uint64_t) 1 << (2 * k)) - 1;
    uint64_t pre_gram = 0;
    uint64_t pre_mask = 0;
    bool canAddNow = false;
    bool isRev = false;

    for (int y = 0; y < len; y++) {

        if (excludeN && (str[y] == 'N' || str[y] == 'n') && (y != lastN + ((k - 1) / 2) - 1)) {
            lastN = y;
        }

        switch (str[y]) {
            case 'A':
            case 'a':
                gram = gram << bits | CODE_A;
                gram_reverse = gram_reverse >> bits | ((uint64_t) CODE_T << ((k - 1) * bits));
                break;
            case 'C':
            case 'c':
                gram = gram << bits | CODE_C;
                gram_reverse = gram_reverse >> bits | ((uint64_t) CODE_G << ((k - 1) * bits));
                break;
            case 'G':
            case 'g':
                gram = gram << bits | CODE_G;
                gram_reverse = gram_reverse >> bits | ((uint64_t) CODE_C << ((k - 1) * bits));
                break;
            case 'T':
            case 't':
                gram = gram << bits | CODE_T;
                gram_reverse = gram_reverse >> bits | ((uint64_t) CODE_A << ((k - 1) * bits));
                break;
            default:
                gram = gram << bits;
                gram_reverse = gram_reverse >> bits;
                break;
        }

        gram = gram & mask;
        gram_reverse = gram_reverse & mask;

        if (!getIdxOnly && lastN <= y - k) {
            if (gram < gram_reverse && (gram & pre_mask) == pre_gram) {
                addingGram = gram;
                isRev = false;
                *getIsRev = false;
                canAddNow = true;
            } else if (gram_reverse < gram && (gram_reverse & pre_mask) == pre_gram) {
                addingGram = gram_reverse;
                isRev = true;
                *getIsRev = true;
                canAddNow = true;
            }
            if (canAddNow) {
                canAddNow = false;
                pos = y - k + 1;
                if ((*ind->kmerToIdx)[addingGram] == 0) { /* A new q-gram */
                    (*ind->kmerToIdx)[addingGram] = kmersCnt + 1;
                    addReadToKmer((*ind->kmerToIdx)[addingGram] - 1, true, rIdx, pos, isRev);
                } else {
                    addReadToKmer((*ind->kmerToIdx)[addingGram] - 1, false, rIdx, pos, isRev);
                }
            }
        }
    }
    if (getIdxOnly) {
        if (gram < gram_reverse && ((gram & pre_mask) == pre_gram) && (*ind->kmerToIdx)[gram] != 0) {
            kIdx = (*ind->kmerToIdx)[gram] - 1;
            *getIsRev = false;
        } else if (gram_reverse < gram && ((gram_reverse & pre_mask) == pre_gram) && (*ind->kmerToIdx)[gram_reverse] != 0) {
            kIdx = (*ind->kmerToIdx)[gram_reverse] - 1;
            *getIsRev = true;
        }
    }
    return kIdx;
}

/**
 * Construct the k-mer index for the reads
 * 
 * @param int      k       The kmer length
 * @param   index
 * @return 1 on success, 0 otherwise
 */
void buildValidKmersHashing(short int k) {
    ind = (struct Indexing *) malloc(sizeof (struct Indexing));
    ind->kmerToIdx = new std::tr1::unordered_map<uint64_t, uint64_t>;

    if (ind->kmerToIdx == NULL) {
        printf("ERROR in allocating a kmerToIdx: Out of memory\n");
        return;
    }
    bool getIsRev;
    for (int i = 0; i < readsCnt; i++) {
        hashKmer(readsPtr[i].seq, readsPtr[i].len, i, true, k, false, &getIsRev);
    }
}

/**
 * A function that frees the allocated memory for kmers
 */
void freeKmers() {
    kmersInvalidCnt = 0;
    for (int i = 0; i < kmersCnt; i++) {
        if (kmersPtr[i].validByKIdx == -1) {
            kmersInvalidCnt += 1;
            if (i < 100) {
                meta << "KIdx = " << i <<"(freq = ("<<kmersPtr[i].rCnt <<","<<kmersPtr[i].revRCnt <<")";
                if (kmersPtr[i].rCnt > 0) {
                    meta << ", RIdx = " << kmersPtr[i].rIds[0]<<", Pos = "<< kmersPtr[i].stPos[0];
                }
                if (kmersPtr[i].revRCnt > 0) {
                    meta << ", RevRIdx = " << kmersPtr[i].revRIds[0]<<", Pos = "<< kmersPtr[i].revStPos[0];
                }
                meta << endl;
            }
        }

        if (kmersPtr[i].rCnt > 0) {
            free(kmersPtr[i].rIds);
            free(kmersPtr[i].stPos);
        }
        if (kmersPtr[i].revRCnt > 0) {
            free(kmersPtr[i].revRIds);
            free(kmersPtr[i].revStPos);
        }
    }
    meta << "kmersInvalidCnt = " << kmersInvalidCnt << endl;
    free(kmersPtr);
    kmersCnt = 0;
    (*ind->kmerToIdx).clear();
    delete ind->kmerToIdx;
    free(ind);
}