#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>

#define SHOW_USAGE() {\
	printf( "\nWelcome to the SHort-read Iterative Trimmer (SHI7) v0.89. Usage:\n");\
	printf( "NINJA_SHI7 in_seqs.fastq out_prefix MINLEN QUAL <QUAL_RIGHT> ...\n");\
	printf( "Choose one of the following optional trimming methods:\n");\
	printf( "   [ROLLING X [PREROLL/POSTROLL]] For rolling average over X bases\n" );\
	printf( "      PREROLL allows initial windows < X. POSTROLL adds back X to both sides.\n");\
	printf( "   [FLOOR X] Cut until X consecutive bases have minimum quality\n");\
	printf( "   [CUT] Treat 'QUAL' and 'QUAL_RIGHT' as number of bases to blindly cut\n");\
	printf( "[ASS_QUALITY X] Chuck reads under minimum average quality X (after trim).\n");\
	printf( "[N-IFY X] Magically transform crap quality (<= X) into ambiguous bases.\n");\
	printf( "[R2 in_seqsR2.fastq]: Invoke paired-end mode with specified R2.\n");\
	printf("\nOh, you'll want to add 31 to all quality args if you're on ol' PHRED64.\n");\
	exit(1); \
}
char BUFFER[4096] = {0};
int main(int argc, char *argv[]) {
	// It's best to start all encounters with arguments.
	if (argc < 5 || argc > 14) SHOW_USAGE()
	int minLen = atoi(argv[3]), tleft = atoi(argv[4]), tright = tleft, doRoll = 0, doCut = 0;
	int preroll = 0, postroll = 0, doFloor = 0, ass_quality = 0, nify = 0; double damned = 0;
	char *pair = 0;
	if (!strcmp(argv[argc-2],"R2"))
		pair = argv[argc-1], argc -= 2,
		printf("Crikey, so these are mate pairs mate.\n");
	if (!strcmp(argv[argc-2],"N-IFY"))
		nify = atoi(argv[argc-1]), argc -=2, 
		printf("N-ifying all bases of quality <= %d.\n",nify);
	if (!strcmp(argv[argc-2],"ASS_QUALITY"))
		ass_quality = atoi(argv[argc-1]), argc -=2, 
		printf("Chucking reads of ass quality (%d).\n",ass_quality);
	if (!strcmp(argv[argc-1],"PREROLL")) 
		--argc, preroll = 1;
	else if (!strcmp(argv[argc-1],"POSTROLL")) 
		--argc, postroll = 1;
	if (!strcmp(argv[argc-2],"ROLLING")) 
		doRoll = atoi(argv[argc-1]), argc -= 2, 
		printf("Using rolling average of %d bases%s.\n",
			doRoll, preroll? " (preroll)" : postroll? " (postroll)" : "");
	else if (!strcmp(argv[argc-2],"FLOOR")) 
		doFloor = atoi(argv[argc-1]), argc -=2,
		printf("Using FLOOR of %d bases.\n",doFloor);
	else if (!strcmp(argv[argc-1],"CUT"))
		printf("Cutting off %d bases on the left and %d on the right.\n",tleft,tright),
		doCut = 1, --argc;
	if (argc > 5)  // Well, the last thing left is QUAL_RIGHT, RIGHT?
		tright = atoi(argv[5]);
	printf("The minimum length a read can be is %d bases.\n",minLen);
	if (!doCut) printf("Trimming while quality < %d on left, < %d on right.\n",tleft,tright);
	
	// Go through their files, C the contents. C what I did there? C?
	sprintf(BUFFER,"%s%s.fq",argv[2],pair? ".R1" : "");
	FILE *in = fopen(argv[1],"rb"), *in2 = 0, *out = fopen(BUFFER,"wb"), *out2 = 0;
	if (pair) 
		sprintf(BUFFER,"%s.R2.fq",argv[2]),
		in2=fopen(pair,"rb"), out2 = fopen(BUFFER,"wb");
	if (!in || !out || (pair && (!in2 || !out2))) {
		puts("Can't open input or output file(s)!"); exit(1);}
	char *Head = malloc(UINT16_MAX), *Seq = malloc(UINT16_MAX), 
		*Shi7 = malloc(UINT16_MAX), *Qual = malloc(UINT16_MAX),
		*Head2 = 0, *Seq2 = 0, *Shi72 = 0, *Qual2 = 0; 
	if (pair) Head2 = malloc(UINT16_MAX), Seq2 = malloc(UINT16_MAX),
		Shi72 = malloc(UINT16_MAX), Qual2 = malloc(UINT16_MAX);
	size_t crap = 0, decency = 0, totI = 0, totJ = 0, totLen = 0, length, length2=0;
	int *Scores = calloc(doRoll,sizeof(*Scores));
	while (Head = fgets(Head,UINT16_MAX,in)) { 
		Seq = fgets(Seq, UINT16_MAX, in);
		Shi7 = fgets(Shi7, UINT16_MAX, in);
		Qual = fgets(Qual, UINT16_MAX, in);
		length = strlen(Qual);
		if (Qual[length-1]=='\n') --length;
		if (Qual[length-1]=='\r') --length;
		if (pair) {
			Head2 = fgets(Head2, UINT16_MAX, in2);
			Seq2 = fgets(Seq2, UINT16_MAX, in2);
			Shi72 = fgets(Shi72, UINT16_MAX, in2);
			Qual2 = fgets(Qual2, UINT16_MAX, in2);
			length2 = strlen(Qual2);
			if (Qual2[length2-1]=='\n') --length2;
			if (Qual2[length2-1]=='\r') --length2;
		}
		int L, R, L2=0, R2=0;
		double darned = 0; // A fitting euphamism
		// Double the mates, half the fun? Quite the mechanism, aPairEndly
		for (int z = 0; z <= (pair != 0); ++z) {
			char *TQual; size_t TLen;
			if (z) TQual = Qual2, TLen = length2;
			else TQual = Qual, TLen = length;
			// Read--no, FEEL the quality (some have gotta 'lookahead').
			int i, j;
			if (doRoll) { // Do the rolling average thing. Fill-all FILO queues
				int total = 0; int roll = 0, pos = 0;
				for (i = 0; i < TLen; i++) {
					Scores[pos] = TQual[i] - 33;
					total += Scores[pos];
					if (++pos == doRoll) pos = 0;
					if (roll == doRoll) total -= Scores[pos];
					else {++roll; if (!preroll) continue;}
					if ((float)total/roll > (float)tleft) break;
				}
				if (postroll) i -= roll - 1;
				if (i == TLen) {++crap; goto HELL;}
				total = 0; roll = 0; pos = 0;
				for (j = TLen-1; j > i; j--) {
					Scores[pos] = TQual[j] - 33;
					total += Scores[pos];
					if (++pos == doRoll) pos = 0;
					if (roll == doRoll) total -= Scores[pos];
					else {++roll; if (!preroll) continue;}
					if ((float)total/roll > (float)tright) break;
				}
				if (postroll) j += roll - 1;
			}
			else if (doFloor) { // Yeah, I'm floored by this sucker too.
				int row = 0;
				for (i = 0; i < TLen; i++) {
					if (TQual[i] - 33 > tleft) ++row;
					else row = 0;
					if (row == doFloor) break;
				}
				i -= doFloor - 1;
				if (i == TLen) {++crap; goto HELL;}
				row = 0;
				for (j = TLen-1; j > i; j--) {
					if (TQual[j] - 33 > tright) ++row;
					else row = 0;
					if (row == doFloor) break;
				}
				j += doFloor - 1;
			} 
			else if (doCut)  // Idiotically slice off a set number of bases.
				i = tright, j = TLen - tleft - 1;
			else { // By the way, the default trimming method is 'go till ya hit somethin good'
				for (i = 0; i < TLen; i++)
					if (TQual[i] - 33 > tleft) break;
				if (i == TLen) {++crap; goto HELL;}
				for (j = TLen-1; j > i; j--)
					if (TQual[j] - 33 > tright) break;
			}
			if (j == i || j - i < minLen) {++crap; goto HELL;}
			++j;
			size_t abs = 0; for (int z = i; z < j; z++) abs += TQual[z];
			double thisAvgQ = (double)abs/(j - i);
			if (thisAvgQ - 33 < ass_quality) {++crap; goto HELL;}
			darned += thisAvgQ;
			if (z) L2 = i, R2 = j;
			else L = i, R = j;
		}
		damned += darned; 
		
		// Regurgitate the sequence(s). They're (prolly) OK enough if they made it this far down.
		fprintf(out,"%s",Head);
		for (int z = L; z < R; z++) fprintf(out,"%c",Qual[z] - 33 > nify ? Seq[z] : 'N');
		fprintf(out,"\n%s",Shi7);
		for (int z = L; z < R; z++) fprintf(out,"%c",Qual[z]);
		fprintf(out,"\n");
		if (pair) {
			fprintf(out2,"%s",Head2);
			for (int z = L2; z < R2; z++) fprintf(out2,"%c",Qual2[z] - 33 > nify ? Seq2[z] : 'N');
			fprintf(out2,"\n%s",Shi72);
			for (int z = L2; z < R2; z++) fprintf(out2,"%c",Qual2[z]);
			fprintf(out2,"\n");
		}
		++decency;
		totLen += (R - L) + (R2 - L2);
		totI += L+L2, totJ += (length - R) + (length2 - R2);

		///// YE SHALL FIND HEREIN THE PITS OF SATAN
		HELL:NULL;
	}
	printf("Done. %ld sequences were utter shi7, quality-wise.\n",crap);
	printf("And %ld sequences were decent.\n", decency);
	printf("The average trimmed length was %0.2f.\n", (double)totLen/decency/(1+(pair!=0)));
	printf("On average, the cut bases were: %0.2f (left), %0.2f (right).\n", 
		(double)totI/decency/(1+(pair!=0)), (double)totJ/decency/(1+(pair!=0)));
	printf("The average read quality was %0.3f.\n-----\n", damned/decency/(1+(pair!=0)) - 33);
	return 0;
}
