
void correctReadsByKmerRead(int *toBeCorrectedRIds, short int *toBeCorrectedStPos, int toBeCorrectedCnt, int correctiveRId, short int correctivePos) {
    int idxInRead = -1;
    uint64_t correctiveFreq = (kmersPtr[kIdxOfMaxFreq].rCnt + kmersPtr[kIdxOfMaxFreq].revRCnt);
    char crrNuc;
    for (int i = 0; i < toBeCorrectedCnt; i++) {
        for (int j = 0; j < nonRelativeLen; j++) {
            idxInRead = toBeCorrectedStPos[i] + nonRelativeVarPos[j];
            if (readsPtr[toBeCorrectedRIds[i]].crrFreq[idxInRead] < correctiveFreq) {
                readsPtr[toBeCorrectedRIds[i]].crrFreq[idxInRead] = correctiveFreq;
                crrNuc = readsPtr[correctiveRId].seq[correctivePos + nonRelativeVarPos[j]];
                if(checkDel){
                    if ((idxInRead - 1) >= 0 && readsPtr[toBeCorrectedRIds[i]].crr[idxInRead - 1] == crrNuc) { //nuc neigh = correct
                        readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = CRR_DEL; //remove
                    } else if ((correctivePos + nonRelativeVarPos[j] - 1) >= 0 
                            && readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] == readsPtr[correctiveRId].seq[correctivePos + nonRelativeVarPos[j] - 1]) { //nuc = correct neigh
                        readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = (char)((int)crrNuc + 100); //insert
                    }  else if ((correctivePos + nonRelativeVarPos[j] + 1) < readsPtr[correctiveRId].len
                            && readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] == readsPtr[correctiveRId].seq[correctivePos + nonRelativeVarPos[j] + 1]) { //nuc = correct neigh
                        readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = (char)((int)crrNuc + 100); //insert
                    } else {
                        readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = crrNuc; //subs
                    }
                }else{
                    readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = crrNuc; //subs
                }
            }
        }
    }
}

void investigateAndCorrectKmer(int toBeCorrectedKIdx) {
    char nuc01 = 'N', revNuc01 = 'N', nuc02 = 'N', revNuc02 = 'N', nuc11 = 'N', revNuc11 = 'N', nuc22 = 'N', revNuc22 = 'N';
    short int checkPos1 = nonRelativeVarPos[0] + 1;
    short int checkPos2 = nonRelativeVarPos[nonRelativeLen - 1] - 1;
    if (kmersPtr[kIdxOfMaxFreq].rCnt > 0) {
        nuc01 = readsPtr[kmersPtr[kIdxOfMaxFreq].rIds[0]].seq[kmersPtr[kIdxOfMaxFreq].stPos[0] + checkPos1];
        nuc11 = readsPtr[kmersPtr[kIdxOfMaxFreq].rIds[0]].seq[kmersPtr[kIdxOfMaxFreq].stPos[0] + checkPos2];
    }
    if (kmersPtr[kIdxOfMaxFreq].revRCnt > 0) {
        revNuc01 = readsPtr[kmersPtr[kIdxOfMaxFreq].revRIds[0]].seq[kmersPtr[kIdxOfMaxFreq].revStPos[0] + checkPos1];
        revNuc11 = readsPtr[kmersPtr[kIdxOfMaxFreq].revRIds[0]].seq[kmersPtr[kIdxOfMaxFreq].revStPos[0] + checkPos2];
    }
    if (kmersPtr[toBeCorrectedKIdx].rCnt > 0) {
        nuc02 = readsPtr[kmersPtr[toBeCorrectedKIdx].rIds[0]].seq[kmersPtr[toBeCorrectedKIdx].stPos[0] + checkPos1];
        nuc22 = readsPtr[kmersPtr[toBeCorrectedKIdx].rIds[0]].seq[kmersPtr[toBeCorrectedKIdx].stPos[0] + checkPos2];
    }
    if (kmersPtr[toBeCorrectedKIdx].revRCnt > 0) {
        revNuc02 = readsPtr[kmersPtr[toBeCorrectedKIdx].revRIds[0]].seq[kmersPtr[toBeCorrectedKIdx].revStPos[0] + checkPos1];
        revNuc22 = readsPtr[kmersPtr[toBeCorrectedKIdx].revRIds[0]].seq[kmersPtr[toBeCorrectedKIdx].revStPos[0] + checkPos2];
    }

    if (nuc01 != 'N' && nuc11 != 'N') {
        if (nuc01 == nuc02 && nuc11 == nuc22) {
            correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].rIds, kmersPtr[toBeCorrectedKIdx].stPos, kmersPtr[toBeCorrectedKIdx].rCnt, kmersPtr[kIdxOfMaxFreq].rIds[0], kmersPtr[kIdxOfMaxFreq].stPos[0]);
        } else if (nuc01 == revNuc02 && nuc11 == revNuc22) {
            correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].revRIds, kmersPtr[toBeCorrectedKIdx].revStPos, kmersPtr[toBeCorrectedKIdx].revRCnt, kmersPtr[kIdxOfMaxFreq].rIds[0], kmersPtr[kIdxOfMaxFreq].stPos[0]);
        }
    }

    if (revNuc01 != 'N' && revNuc11 != 'N') {
        if (revNuc01 == nuc02 && revNuc11 == nuc22) {
            correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].rIds, kmersPtr[toBeCorrectedKIdx].stPos, kmersPtr[toBeCorrectedKIdx].rCnt, kmersPtr[kIdxOfMaxFreq].revRIds[0], kmersPtr[kIdxOfMaxFreq].revStPos[0]);
        } else if (revNuc01 == revNuc02 && revNuc11 == revNuc22) {
            correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].revRIds, kmersPtr[toBeCorrectedKIdx].revStPos, kmersPtr[toBeCorrectedKIdx].revRCnt, kmersPtr[kIdxOfMaxFreq].revRIds[0], kmersPtr[kIdxOfMaxFreq].revStPos[0]);
        }
    }
}

void correctKmerByMaxFreqKmer(int toBeCorrectedKIdx, bool toBeCorrectedKIsRev, int kIdxOfMaxFreq, bool kOfMaxFreqIsRev) {
    bool notCorrectedYet = false;
    if (!kOfMaxFreqIsRev) {
        if (kmersPtr[kIdxOfMaxFreq].rCnt > 0) {
            if (!toBeCorrectedKIsRev) {
                if (kmersPtr[toBeCorrectedKIdx].rCnt > 0) {
                    correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].rIds, kmersPtr[toBeCorrectedKIdx].stPos, kmersPtr[toBeCorrectedKIdx].rCnt, kmersPtr[kIdxOfMaxFreq].rIds[0], kmersPtr[kIdxOfMaxFreq].stPos[0]);
                } else {
                    notCorrectedYet = true;
                }
            } else {
                if (kmersPtr[toBeCorrectedKIdx].revRCnt > 0) {
                    correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].revRIds, kmersPtr[toBeCorrectedKIdx].revStPos, kmersPtr[toBeCorrectedKIdx].revRCnt, kmersPtr[kIdxOfMaxFreq].rIds[0], kmersPtr[kIdxOfMaxFreq].stPos[0]);
                } else {
                    notCorrectedYet = true;
                }
            }
        } else {
            notCorrectedYet = true;
        }
    } else {

        if (kmersPtr[kIdxOfMaxFreq].revRCnt > 0) {
            if (!toBeCorrectedKIsRev) {
                if (kmersPtr[toBeCorrectedKIdx].rCnt > 0) {
                    correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].rIds, kmersPtr[toBeCorrectedKIdx].stPos, kmersPtr[toBeCorrectedKIdx].rCnt, kmersPtr[kIdxOfMaxFreq].revRIds[0], kmersPtr[kIdxOfMaxFreq].revStPos[0]);
                } else {
                    notCorrectedYet = true;
                }
            } else {
                if (kmersPtr[toBeCorrectedKIdx].revRCnt > 0) {
                    correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].revRIds, kmersPtr[toBeCorrectedKIdx].revStPos, kmersPtr[toBeCorrectedKIdx].revRCnt, kmersPtr[kIdxOfMaxFreq].revRIds[0], kmersPtr[kIdxOfMaxFreq].revStPos[0]);
                } else {
                    notCorrectedYet = true;
                }
            }
        } else {
            notCorrectedYet = true;
        }
    }
    if (notCorrectedYet) {
        investigateAndCorrectKmer(toBeCorrectedKIdx);
    }
}

void correctKmerGps() {
    int i, j;
    for (i = 0; i < kmerGpsCnt; i++) {
        for(j = 0; j < kmerGpsPtr[i].kCnt; j++){
            if (j != kmerGpsPtr[i].bestKIdx){
                correctKmerByMaxFreqKmer(kmerGpsPtr[i].kIds[j], kmerGpsPtr[i].isRev[j], kmerGpsPtr[i].kIds[kmerGpsPtr[i].bestKIdx], kmerGpsPtr[i].isRev[kmerGpsPtr[i].bestKIdx]);
            }
        }
    }
}

void updateReadsByCorrections() {
    short int oldLen = 0;
    string oldSeq, oldQV;
    short int y;
    for (int i = 0; i < readsCnt; i++) {
        oldLen = readsPtr[i].len;
        oldSeq = readsPtr[i].seq;
        oldQV = readsPtr[i].qv;
        //Setting the corrected sequence
        readsPtr[i].len = 0;
        for (y = 0; y < oldLen; y++) {
            if (readsPtr[i].crr[y] != CRR_DEL) {
                readsPtr[i].seq = (char *) realloc(readsPtr[i].seq, (readsPtr[i].len + 1) * sizeof (char));
                readsPtr[i].qv = (char *) realloc(readsPtr[i].qv, (readsPtr[i].len + 1) * sizeof (char));
                if (readsPtr[i].crr[y] == CRR_NO) {
                    readsPtr[i].seq[readsPtr[i].len] = oldSeq[y];
                    readsPtr[i].qv[readsPtr[i].len] = oldQV[y];
                } else if (readsPtr[i].crr[y] >= CRR_SUBS_A && readsPtr[i].crr[y] <= CRR_SUBS_t) {
                    readsPtr[i].seq[readsPtr[i].len] = readsPtr[i].crr[y];
                    readsPtr[i].qv[readsPtr[i].len] = (oldQV[y] < 'c') ? 'c' : oldQV[y];
                } else {// if ((int)(readsPtr[i].crr[y] - 100) >= CRR_SUBS_A && (int)(readsPtr[i].crr[y] - 100) <= CRR_SUBS_t) {
                    readsPtr[i].seq[readsPtr[i].len] = (char) ((int) readsPtr[i].crr[y] - 100);
                    readsPtr[i].qv[readsPtr[i].len] = (oldQV[y] < 'c') ? 'c' : oldQV[y];
                    readsPtr[i].len += 1;
                    readsPtr[i].seq = (char *) realloc(readsPtr[i].seq, (readsPtr[i].len + 1) * sizeof (char));
                    readsPtr[i].qv = (char *) realloc(readsPtr[i].qv, (readsPtr[i].len + 1) * sizeof (char));
                    readsPtr[i].seq[readsPtr[i].len] = oldSeq[y];
                    readsPtr[i].qv[readsPtr[i].len] = oldQV[y];
                }
                readsPtr[i].len += 1;
            }
        }
        readsPtr[i].seq = (char *) realloc(readsPtr[i].seq, (readsPtr[i].len + 1) * sizeof (char));
        readsPtr[i].qv = (char *) realloc(readsPtr[i].qv, (readsPtr[i].len + 1) * sizeof (char));
        readsPtr[i].seq[readsPtr[i].len] = '\0';
        readsPtr[i].qv[readsPtr[i].len] = '\0';
        //Resetting corrections
        readsPtr[i].crr = (char *) realloc(readsPtr[i].crr, (readsPtr[i].len) * sizeof (char));
        readsPtr[i].crrFreq = (int *) realloc(readsPtr[i].crrFreq, (readsPtr[i].len) * sizeof (int));
        if (readsPtr[i].crr == 0 || readsPtr[i].crrFreq == 0) {
            printf("ERROR in reallocating a readsPtr[i].crr/readsPtr[i].crrFreq: Out of memory\n");
            return;
        }
        for (y = 0; y < readsPtr[i].len; y++) {
            readsPtr[i].crr[y] = CRR_NO;
            readsPtr[i].crrFreq[y] = 0;
        }
    }
}
/*
TP = 16,902,900
FP = 354,564
FN = 9,907,200
TN = 130,157,136
Sensitivity = 63.0468%
Specificity = 99.7283%
Accuracy = 61.7243%
Gain = 61.7243%
        (16,902,900+130,157,136)/(16,902,900+130,157,136+9,907,200+354,564) = 0.93477214219
        
TP = 16,904,268
FP = 355,536
FN = 9,905,832
TN = 130,156,164
Sensitivity = 63.0519%
Specificity = 99.7276%
Accuracy = 61.7257%
Gain = 61.7257%
                (16,904,268+130,156,164)/(16,904,268+130,156,164+9,905,832+355,536) = 0.93477465932	
        
TP = 13,670,784
FP = 542,880
FN = 13,139,316
TN = 129,968,820
Sensitivity = 50.9912%
Specificity = 99.584%
Accuracy = 48.9663%
Gain = 48.9663%
(13,670,784+129,968,820)/(13,670,784+13,139,316+542,880+129,968,820)    = 0.91303051452	    
*/
