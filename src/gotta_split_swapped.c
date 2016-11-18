#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#define LINELEN 1000000
#if WIN32 || _WIN32 || __WIN32__
#include <windows.h>
#endif
#include <time.h>
#define PRINT_USAGE() { \
	printf( "\nG0TTA-SPLI7 v0.93s by Gabe. Splits PE FASTQ by MOTHUR oligos.txt\n"); \
	printf( "Usage: gotta_split r1.fq r2.fq oligos.txt [PREFIX] [--deep] [--errors N]\n");\
	printf( "Add a PREFIX to append to the filename. Can be \"dirname/\"\n"); \
	printf( "--deep also looks inside the reads, not just the beginning.\n"); \
	printf( "(Only supports BARCODE-exclusive oligos.txt w/only N's as ambigs)\n"); \
	exit(1);\
}
uint8_t *EXPAND[128] = {0}, *MUTATE[128] = {0};
typedef struct SampObj SampObj;
struct SampObj { 
	SampObj *fLink, *rLink;
	FILE *fileF, *fileR;
	char *tagF, *tagR, *name;
}; 
static int descSampCmpF(const void *restrict a, const void *restrict b) {
	SampObj *s1 = *(SampObj **)a, *s2 = *(SampObj **)b;
	int fcomp = strcmp(s2->tagF, s1->tagF);
	if (fcomp) return fcomp;
	//return s2 < s1 ? -1 : s2 > s1;
	return (int64_t)s2 - (int64_t)s1;
}
static int descSampCmpFL(const void *restrict a, const void *restrict b) {
	SampObj *s1 = *(SampObj **)a, *s2 = *(SampObj **)b;
	return strcmp(s2->tagF, s1->tagF);
}
static int descSampCmpR(const void *restrict a, const void *restrict b) {
	SampObj *s1 = *(SampObj **)a, *s2 = *(SampObj **)b;
	int fcomp = strcmp(s2->tagR, s1->tagR);
	if (fcomp) return fcomp;
	//return s2 < s1 ? -1 : s2 > s1;
	return (int64_t)s2 - (int64_t)s1;
}
static int descSampCmpRL(const void *restrict a, const void *restrict b) {
	SampObj *s1 = *(SampObj **)a, *s2 = *(SampObj **)b;
	return strcmp(s2->tagR, s1->tagR);
}

void addDictIn(char ***Dicts, SampObj ***SCaches, size_t *Dix, size_t *Dsz, 
 SampObj *samp, char *word, int wix, int ne, int maxE, char *newEntry) {
	if (!word[wix]) {
		if (Dix[ne] == Dsz[ne]) { // resize dict if necessary
			Dicts[ne] = realloc(Dicts[ne],sizeof(*Dicts[ne])*(Dsz[ne]*=2));
			SCaches[ne] = realloc(SCaches[ne],sizeof(*SCaches[ne])*Dsz[ne]);
			if (!Dicts[ne] || !SCaches[ne]) {puts("Out of dictionary memory."); exit(3);}
		}
		newEntry[wix] = 0;
		Dicts[ne][Dix[ne]] = realloc(newEntry, wix+1);
		SCaches[ne][Dix[ne]] = samp;
		++Dix[ne];
		return;
	}
	if (ne < maxE && MUTATE[word[wix]]) { // do errors
		int toMutate = *MUTATE[word[wix]] - 47;
		for (int i = 1; i < toMutate; ++i) {
			char *newWrd = malloc(1024);
			for (int j = 0; j < wix; ++j) newWrd[j] = newEntry[j];
			newWrd[wix] = MUTATE[word[wix]][i];
			addDictIn(Dicts, SCaches, Dix, Dsz, samp, word, wix+1, ne+1, maxE, newWrd);
		}
	}
	if (EXPAND[word[wix]]) { // if not base letter, expand all
		int toExpand = *EXPAND[word[wix]] - 47;
		for (int i = 1; i < toExpand; ++i) {
			char *newWrd = malloc(1024);
			for (int j = 0; j < wix; ++j) newWrd[j] = newEntry[j];
			newWrd[wix] = EXPAND[word[wix]][i];
			addDictIn(Dicts, SCaches, Dix, Dsz, samp, word, wix+1, ne, maxE, newWrd);
		}
	} 
	else newEntry[wix] = word[wix], addDictIn(Dicts, SCaches, Dix, Dsz, samp, word, 
		wix+1, ne, maxE, newEntry);
}
static int strPcmp(const void *a, const void *b) {
	return strcmp(**(char ***)a,**(char ***)b); }

// Technically you can make it return char *p, on which you can manually do p - String
inline size_t crBSTw(char *key, size_t sz, char **String, char **what) {
	char **p = String;
	while (sz) {
		size_t w = sz >> 1;
		char *ref_s = *(p+w+1), *key_s = key;
		while (*ref_s == *key_s && *ref_s && *key_s) ++ref_s, ++key_s; 
		if (!*ref_s) { *what = key_s; return p+w+1-String; }
		if (*ref_s < *key_s) p+=w+1, sz-=w+1;
		//while (*ref_s == *key_s++) if (!*ref_s++) { *what = key_s-1; return p+w+1-String; }
		//if (*ref_s < *(key_s-1)) { p+=w+1; sz-=w+1; }
		else sz = w;
	}
	char *ref_s = *p, *key_s = key;
	//while (*ref_s == *key_s && *ref_s && *key_s) ++ref_s, ++key_s; 
	//if (!*ref_s) { *what = key_s; return p - String; }
	while (*ref_s++ == *key_s) if (!*key_s++) return -1;
	if (!*(ref_s-1)) { *what = key_s; return p - String; }
	//*what = 0; 
	return -1;
}
int main( int argc, char *argv[] )
{
	clock_t start = clock();
	if ( argc < 4 || argc > 8 ) PRINT_USAGE()
	#if WIN32 || _WIN32 || __WIN32__
	_setmaxstdio(2048); // this doesn't work on *Nix...
	#endif
	char fname[4096], *prefix = "", doInread = 0;
	int errors = 0;
	if (argc > 5 && !strcmp(argv[argc-2],"--errors")) 
		--argc, errors = atoi(argv[argc--]), printf("Allowing %d errors\n",errors);
	if (argc > 4 && !strcmp(argv[argc-1],"--deep"))
		--argc, doInread = 1, puts("Looking within reads...");
	if (argc > 4)
		--argc, prefix = argv[4], printf("Using prefix: %s\n",prefix);
	FILE *r1 = fopen(argv[1], "rb"), *r2 = fopen(argv[2], "rb"), *olig = fopen(argv[3], "rb");
	if (!r1 || !r2 || !olig) { fputs("Can't open input(s)!\n",stderr); exit(1); }
	
	//Parse the oligos.txt file
	fseek(olig, 0, SEEK_END); size_t sz = ftell(olig); rewind(olig);
	char *OligoDump = malloc(sz+1);
	fread(OligoDump,1,sz,olig); OligoDump[sz] = 0;
	int numSamps = 0; for (int i = 0; i < sz; ++i) numSamps += OligoDump[i]=='\n';
	if (OligoDump[sz-1] != '\n') ++numSamps;
	printf("Number of samples = %d\n",numSamps);
	if (numSamps > 1000) puts("Warning. Your system may not handle over 2000 files");
	SampObj **Samps = malloc(numSamps*sizeof(*Samps));
	char *OligoPtr = OligoDump;
	/* for (int i = 0; i < numSamps; ++i) {
		Samps[i] = calloc(1,sizeof(*Samps[i]));
		while (*OligoPtr++ != '\t'); Samps[i]->tagF = OligoPtr;
		while (*++OligoPtr != '\t'); *OligoPtr++ = 0, Samps[i]->tagR = OligoPtr;
		while (*++OligoPtr != '\t'); *OligoPtr++ = 0, Samps[i]->name = OligoPtr; 
		while (*OligoPtr != '\n' && *OligoPtr != '\r') ++OligoPtr; *OligoPtr = 0;
		sprintf(fname,"%s%s_R1.fastq",prefix,Samps[i]->name);
		Samps[i]->fileF = fopen(fname,"wb");;
		if (!Samps[i]->fileF) { printf("can't open %s\n",fname); exit(1); }
		sprintf(fname,"%s%s_R2.fastq",prefix,Samps[i]->name);
		Samps[i]->fileR = fopen(fname,"wb");
		if (!Samps[i]->fileR) { printf("can't open %s\n",fname); exit(1); }
	} */
	rewind(olig);
	char *sampLine = malloc(LINELEN+1), *sampLineO = sampLine;
	for (int i = 0; i < numSamps; ++i) {
		sampLine = fgets(sampLine,LINELEN,olig);
		if (!sampLine) {printf("Critical error on oligo read-in, line %i\n",i); exit(2);}
		int j=0; 
		for (; j < LINELEN; j++) {
			if (!sampLine[j] || sampLine[j] == '\n' || sampLine[j] == '\r') 
				{printf("ERROR. Line %d parses wrongly.\n",i); exit(2);}
			if (sampLine[j] == '\t') break;
		}
		if (j == LINELEN) {printf("ERROR. Post-line %d parses wrongly.\n",i); exit(2);}
		int start = ++j;
		for (; j < LINELEN; j++) {
			if (!sampLine[j] || sampLine[j] == '\n' || sampLine[j] == '\r') 
				{printf("ERROR. Line %d parses wrongly.\n",i); exit(2);}
			if (sampLine[j] == '\t') break;
		}
		if (j == LINELEN) {printf("ERROR. Post-line %d parses wrongly.\n",i); exit(2);}
		char *field = malloc(j - start + 1);
		memcpy(field, sampLine + start, j - start);
		field[j-start] = 0;
		//if (j-start < minTag) minTag = j-start; // figure out if this length is smaller than expected
		
		start = ++j;
		for (; j < LINELEN; j++) {
			if (!sampLine[j] || sampLine[j] == '\n' || sampLine[j] == '\r') 
				{printf("ERROR. Line %d parses wrongly.\n",i); exit(2);}
			if (sampLine[j] == '\t') break;
		}
		if (j == LINELEN) {printf("ERROR. Post-line %d parses wrongly.\n",i); exit(2);}
		char *field2 = malloc(j - start + 1);
		memcpy(field2, sampLine + start, j - start);
		field2[j-start] = 0;
		//if (j-start < minTag) minTag = j-start; 
		
		start = ++j;
		for (; j < LINELEN; j++) {
			if (!sampLine[j] || sampLine[j] == '\n' || sampLine[j] == '\r') 
				break; // ok to see these tokens after this field
			if (sampLine[j] == '\t') break;
		}
		if (j == LINELEN) {printf("ERROR. Post-line %d parses wrongly.\n",i); exit(2);}
		char *field3 = malloc(j - start + 1);
		memcpy(field3, sampLine + start, j-start);
		field3[j-start] = 0;
		
		sprintf(fname,"%s%s_R1.fastq\0",prefix,field3);
		FILE *temp = fopen(fname,"wb"); //Samps[i]->fileF = temp;
		if (!temp) { printf("can't open %s\n",fname); exit(1); }
		sprintf(fname,"%s%s_R2.fastq\0",prefix,field3);
		FILE *temp2 = fopen(fname,"wb"); //Samps[i]->fileR = temp2;
		if (!temp2) { printf("can't open %s\n",fname); exit(1); }
		
		Samps[i] = malloc(sizeof(*Samps[i])); 
		/* struct SampObj { 
			SampObj *fLink, *rLink;
			FILE *fileF, *fileR;
			char *tagF, *tagR, *name;
		};  */
		*Samps[i] = (SampObj) {0,0, temp, temp2, field, field2, field3};
		//printf("Sample %d: %s, %s, %s\n",i, Samps[i]->name, Samps[i]->tagF, Samps[i]->tagR);
	}
	printf("Parsed file.\n");
	//printf("test complete. Stopping...\n"); exit(1);
	//populate the nucleotide redundancy expander. Currently supports only N, but trivial to add more.
	for (int i = 0; i < 128; ++i) EXPAND[i] = 0; // set null for no expansion of this letter
	EXPAND['N'] = "4ACGT\0"; // EXPAND['R']= "3ACG\0", etc
	EXPAND['n'] = EXPAND['N'];
	MUTATE['A'] = "3CGT\0", MUTATE['C'] = "3AGT\0", MUTATE['G'] = "3ACT\0", 
		MUTATE['T'] = "3ACG\0";
	MUTATE['a'] = MUTATE['A'], MUTATE['c'] = MUTATE['C'], MUTATE['g'] = MUTATE['G'],
		MUTATE['t'] = MUTATE['T'];
	
	// Set up dictionary (forward and reverse)
	size_t *DictSzF = malloc((errors+1)*sizeof(*DictSzF)), 
		*DictSzR = malloc((errors+1)*sizeof(*DictSzR)),
		*DixF = calloc(errors+1,sizeof(*DixF)),
		*DixR = calloc(errors+1,sizeof(*DixR));
	char ***DictF = malloc((errors+1)*sizeof(*DictF)),
		***DictR = malloc((errors+1)*sizeof(*DictR));
	SampObj ***SCacheF = malloc((errors+1)*sizeof(*SCacheF)),
		***SCacheR = malloc((errors+1)*sizeof(*SCacheR));
	int initSz = 100;
	for (int i = 0; i <= errors; ++i) {
		DictSzF[i] = initSz, DictSzR[i] = initSz,
		DictF[i] = malloc(initSz*sizeof(*DictF[i])),
		DictR[i] = malloc(initSz*sizeof(*DictR[i])),
		SCacheF[i] = malloc(initSz*sizeof(*SCacheF[i])),
		SCacheR[i] = malloc(initSz*sizeof(*SCacheR[i]));
	}
	
	//addDictIn(char ***Dicts, SampObj ***SCaches, size_t *Dix, size_t *Dsz, 
	//   SampObj *samp, char *word, int wix, int ne, int maxE, char *newEntry)
	// deduplicate Samps by chaining to one another, expanding and adding to dict. 
	qsort(Samps,numSamps,sizeof(*Samps),descSampCmpF);
	for (int i = 1; i < numSamps; ++i) {
		if (!descSampCmpFL(Samps+i-1,Samps+i)) // same, so chain them
			Samps[i]->fLink = Samps[i-1];
		else // expand, dedupe, etc on last 
			addDictIn(DictF, SCacheF, DixF, DictSzF, Samps[i-1], Samps[i-1]->tagF, 0, 0, 
				errors, malloc(1024));
	}
	addDictIn(DictF, SCacheF, DixF, DictSzF, Samps[numSamps-1], Samps[numSamps-1]->tagF, 
		0, 0, errors, malloc(1024));
	
	qsort(Samps,numSamps,sizeof(*Samps),descSampCmpR);
	for (int i = 1; i < numSamps; ++i) {
		if (!descSampCmpRL(Samps+i-1,Samps+i)) // same, so chain them
			Samps[i]->rLink = Samps[i-1];
		else // expand, dedupe, etc on last 
			addDictIn(DictR, SCacheR, DixR, DictSzR, Samps[i-1], Samps[i-1]->tagR, 0, 0,
				errors, malloc(1024));
	}
	addDictIn(DictR, SCacheR, DixR, DictSzR, Samps[numSamps-1], Samps[numSamps-1]->tagR, 
		0, 0, errors, malloc(1024));
	printf("Number of things added F: %llu [%llu], R: %llu[%llu]\n",DixF[0],DixF[1],DixR[0],DixR[1]);

	// sort the dictionary at all error levels and copy to new array
	char ***DictSrtF = malloc((errors+1)*sizeof(*DictSrtF));
	SampObj ***SOSrtF = malloc((errors+1)*sizeof(*SOSrtF));
	for (int i = 0; i <= errors; ++i) {
		char ***DictPtr = malloc(DixF[i]*sizeof(*DictPtr));
		for (size_t j = 0; j < DixF[i]; ++j) DictPtr[j] = DictF[i] + j;
		qsort(DictPtr,DixF[i],sizeof(*DictPtr),strPcmp);
		DictSrtF[i] = malloc(DixF[i]*sizeof(*DictSrtF[i]));
		SOSrtF[i] = malloc(DixF[i]*sizeof(*SOSrtF[i]));
		for (size_t j = 0; j < DixF[i]; ++j) DictSrtF[i][j] = *DictPtr[j], 
			SOSrtF[i][j] = SCacheF[i][DictPtr[j]-DictF[i]];
		free(DictF[i]); free(SCacheF[i]); free(DictPtr);
	}
	char ***DictSrtR = malloc((errors+1)*sizeof(*DictSrtR));
	SampObj ***SOSrtR = malloc((errors+1)*sizeof(*SOSrtR));
	for (int i = 0; i <= errors; ++i) {
		char ***DictPtr = malloc(DixR[i]*sizeof(*DictPtr));
		for (size_t j = 0; j < DixR[i]; ++j) DictPtr[j] = DictR[i] + j;
		qsort(DictPtr,DixR[i],sizeof(*DictPtr),strPcmp);
		DictSrtR[i] = malloc(DixR[i]*sizeof(*DictSrtR[i]));
		SOSrtR[i] = malloc(DixR[i]*sizeof(*SOSrtR[i]));
		for (size_t j = 0; j < DixR[i]; ++j) DictSrtR[i][j] = *DictPtr[j], 
			SOSrtR[i][j] = SCacheR[i][DictPtr[j]-DictR[i]];
		free(DictR[i]); free(SCacheR[i]); free(DictPtr);
	}
	
	
	
	/* char ***DictPtr = malloc(dictSzF*sizeof(*DictPtr));
	for (size_t i = 0; i < dictSzF; ++i) DictPtr[i] = DictF + i;
	qsort(DictPtr,dictSzF,sizeof(*DictPtr),strPcmp);
	char **DictSrtF = malloc(dictSzF*sizeof(*DictSrtF));
	SampObj **SOSrtF = malloc(dictSzF*sizeof(*SOSrtF));
	for (size_t i = 0; i < dictSzF; ++i) DictSrtF[i] = *DictPtr[i], SOSrtF[i] = SCacheF[DictPtr[i]-DictF];
	free(DictF); free(SCacheF); free(DictPtr); 
	
	
	DictPtr = malloc(dictSzR*sizeof(*DictPtr));
	for (size_t i = 0; i < dictSzR; ++i) DictPtr[i] = DictR + i;
	qsort(DictPtr,dictSzR,sizeof(*DictPtr),strPcmp);
	char **DictSrtR = malloc(dictSzR*sizeof(*DictSrtR));
	SampObj **SOSrtR = malloc(dictSzR*sizeof(*SOSrtR));
	for (size_t i = 0; i < dictSzR; ++i) DictSrtR[i] = *DictPtr[i], SOSrtR[i] = SCacheR[DictPtr[i]-DictR];
	free(DictR); free(SCacheR); free(DictPtr);
	*/
	
	// Go through the fastq's and assign one by one (somebody needs to multithread this)
	char *Head = malloc(LINELEN+1), *Seq = malloc(LINELEN+1), 
		*Shi7 = malloc(LINELEN+1), *Qual = malloc(LINELEN+1),
		*Head2 = malloc(LINELEN+1), *Seq2 = malloc(LINELEN+1), 
		*Shi72 = malloc(LINELEN+1), *Qual2 = malloc(LINELEN+1),
		*whereF, *whereR; 
	size_t numSeqs = 0, badSeqs = 0;
	#ifdef LOG
	FILE *log = fopen("split.log","wb");
	if (!log) puts("Can't open log. Oh well.");
	#endif
	//printf("test complete. Stopping...\n"); exit(1);
	while (Head = fgets(Head,LINELEN,r1)) { 
		++numSeqs;
		Seq = fgets(Seq, LINELEN, r1);
		Shi7 = fgets(Shi7, LINELEN, r1);
		Qual = fgets(Qual, LINELEN, r1);
		
		Head2 = fgets(Head2, LINELEN, r2);
		Seq2 = fgets(Seq2, LINELEN, r2);
		Shi72 = fgets(Shi72, LINELEN, r2);
		Qual2 = fgets(Qual2, LINELEN, r2);
		char *r1, *r2, *q1, *q2;
		
		#ifdef LOG
		fprintf(log,"R1=%sR2=%s",Seq,Seq2);
		#endif
		
		//for (int x = 0; x <= errors; ++x) {
		int x = -1, y = -1; // x is dict/cache, y is left
		size_t ixF = -1, ixR = -1;
		while (ixF == (size_t)-1 && ++x <= errors) {
			r1 = Seq,  q1 = Qual;
			ixF = crBSTw(r1,DixF[x]-1,DictSrtF[x],&whereF);
			if (doInread) 
				while (ixF == (size_t)-1 && *r1) ++r1, ++q1, ixF = crBSTw(r1,DixF[x]-1,DictSrtF[x],&whereF);
		}
		if (ixF != (size_t)-1) while (ixR == (size_t)-1 && ++y <= errors) {
			r2 = Seq2, q2 = Qual2;
			ixR = crBSTw(r2,DixR[y]-1,DictSrtR[y],&whereR);
			if (doInread)
				while (ixR == (size_t)-1 && *r2) ++r2, ++q2, ixR = crBSTw(r2,DixR[y]-1,DictSrtR[y],&whereR);
		}
		
		if (ixF != (size_t)-1 && ixR != (size_t)-1) {
			SampObj *f = SOSrtF[x][ixF], *r = SOSrtR[y][ixR];
			while (f != r) { // LOCK-STEP ALGORITHM
				while (f && f < r) f = f->fLink;
				if (!f) break;
				while (r && r < f) r = r->rLink;
				if (!r) break;
			}
			if (f==r) {
				#ifdef LOG
				fprintf(log,"-->Concordant; Sample = %s [%d,%d]\n",f->name, x,y);
				#endif
				size_t offset = whereF - r1;
				fprintf(f->fileF,"%s%s%s%s",Head, r1+offset, Shi7, q1+offset);
				offset = whereR - r2;
				fprintf(f->fileR,"%s%s%s%s",Head2, r2+offset, Shi72, q2+offset);
				goto COMMIT;
			}
		}
		else {
			x = -1, y = -1, ixF = -1, ixR = -1;
			while (ixF == (size_t)-1 && ++y <= errors) {
				r1 = Seq, q1 = Qual;
				ixF = crBSTw(r1,DixR[y]-1,DictSrtR[y],&whereF);
				if (doInread) 
					while (ixF == (size_t)-1 && *r1) ++r1, ++q1, ixF = crBSTw(r1,DixR[y]-1,DictSrtR[y],&whereF);
			}
			if (ixF != (size_t)-1) while (ixR == (size_t)-1 && ++x <= errors) {
				r2 = Seq2, q2 = Qual2;
				ixR = crBSTw(r2,DixF[x]-1,DictSrtF[x],&whereR);
				if (doInread) 
					while (ixR == (size_t)-1 && *r2) ++r2, ++q2, ixR = crBSTw(r2,DixF[x]-1,DictSrtF[x],&whereR);
			}
			
			if (ixF != (size_t)-1 && ixR != (size_t)-1) {
				SampObj *f = SOSrtR[y][ixF], *r = SOSrtF[x][ixR];
				while (f != r) { // LOCK-STEP ALGORITHM
					while (f && f < r) f = f->rLink;
					if (!f) break;
					while (r && r < f) r = r->fLink;
					if (!r) break;
				}
				if (f==r) {
					#ifdef LOG
					fprintf(log,"-->DISCORDANT; Sample = %s [%d,%d]\n",f->name, x, y);
					#endif
					size_t offset = whereF - r1;
					fprintf(f->fileR,"%s%s%s%s",Head, r1+offset, Shi7, q1+offset);
					offset = whereR - r2;
					fprintf(f->fileF,"%s%s%s%s",Head2, r2+offset, Shi72, q2+offset);
					goto COMMIT;
				}
			}
		}
		//}
		++badSeqs;
		#ifdef LOG
		fputs("NO ALIGNMENT FOUND\n",log);
		#endif
		COMMIT: NULL;
		#ifdef LOG
		fputc('\n',log);
		#endif
	}
	printf("Done. Processed %llu sequences (%llu failures).\n",numSeqs,badSeqs);
	printf("Time: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock(); 
}
