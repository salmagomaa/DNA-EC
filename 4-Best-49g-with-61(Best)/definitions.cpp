//Definitions
#define MAX_NUM_QGRAMS 1048576 /* The maximum number of q-grams. Used for allocating some initial memory structures. */
#define INVALID ((uint64_t)0xffffffffffffffff)
#define DEFAULT_Q 21  /* The default length of a q-gram */
#define CODE_A 0x00
#define CODE_C 0x01
#define CODE_G 0x02
#define CODE_T 0x03
#define CRR_NO      33    //noCrr{!}
#define CRR_DEL     45    //Del{-}
#define CRR_INS_A   165   //(65+100)->Ins{A}
#define CRR_INS_a   197   //(97+100)->Ins{a}
#define CRR_INS_C   167   //(67+100)->Ins{C}
#define CRR_INS_c   199   //(99+100)->Ins{c}
#define CRR_INS_G   171   //(71+100)->Ins{G}
#define CRR_INS_g   203   //(103+100)->Ins{g}
#define CRR_INS_T   184   //(84+100)->Ins{T}
#define CRR_INS_t   216   //(116+100)->Ins{t}
#define CRR_SUBS_A  65    //Subs{A}
#define CRR_SUBS_a  97    //Subs{a}
#define CRR_SUBS_C  67    //Subs{C}
#define CRR_SUBS_c  99    //Subs{c}
#define CRR_SUBS_G  71    //Subs{G}
#define CRR_SUBS_g  103   //Subs{g}
#define CRR_SUBS_T  84    //Subs{T}
#define CRR_SUBS_t  116   //Subs{t}

//Structures
struct Read {
    short int len;
    char *seq;
    char *qv;
    char *crr;
    int *crrFreq;
};

struct Kmer {
    int rCnt;
    int *rIds;
    short int *stPos;
    int revRCnt;
    int *revRIds;
    short int *revStPos;
    uint64_t validByKIdx; //-1 ->Not validated yet, itself->is valid, otherwise-> idx of its valid kmer idx
};

struct Indexing {
    std::tr1::unordered_map<uint64_t, uint64_t> *kmerToIdx; /* Mapping from solid-gram to its id */
};

//struct KmerGp {
//    uint64_t kCnt;
//    uint64_t *kIds;
//    uint64_t bestKIdx;
//};
//
//Global Variables
short int seqLn;
struct Read *readsPtr;
int readsCnt = 0;
struct Kmer *kmersPtr;
uint64_t kmersCnt = 0;
struct Indexing *ind; /* The k-mer index */
int kmersInvalidCnt = 0;
int *corresKIdxs;
char **corresKStrs;
bool *corresKIsRev;
int corresKStrCnt = 0;
int corresKIdxCnt = 0;
int kIdxOfMaxFreq = -1;
bool kOfMaxFreqIsRev = false;
bool checkDel;
//struct KmerGp *kGpsPtr;
//uint64_t kGpsCnt = 0;
//short int k;
//short int loops;
////short int maxUnMatCnt = 3;
//
//
//
