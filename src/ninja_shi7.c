#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

#define SHOW_USAGE() {\
	printf( "\nNINJA Is Not Just Another - SHort-read Iterative Trimmer. Usage:\n");\
	printf( "NINJA_SHI7 in_seqs.fastq out_seqs.fastq MINLEN QUAL <QUAL_RIGHT> ...\n");\
	printf( "Choose one of the following optional trimming methods:\n");\
	printf( "   [ROLLING X [PREROLL/POSTROLL]] For rolling average over X bases\n" );\
	printf( "      PREROLL allows initial windows < X. POSTROLL adds back X to both sides.\n");\
	printf( "   [FLOOR X] Cut until X consecutive bases have minimum quality\n");\
	printf( "   [CUT] Treat 'QUAL' and 'QUAL_RIGHT' as number of bases to blindly cut\n");\
	printf( "[ASS_QUALITY X] Chuck reads under minimum average quality X (post-trim).\n");\
	printf( "[N-IFY X] Magically transform crap quality into ambiguous bases.\n");\
	printf("\nOh, you'll want to add 31 to your qualities if you're on ol' PHRED64.\n");\
	exit(1); \
}

int main(int argc, char *argv[]) {
	// It's best to start all encounters with arguments.
	if (argc < 5 || argc > 14) SHOW_USAGE()
	int minLen = atoi(argv[3]), tleft = atoi(argv[4]), tright = tleft, doRoll = 0, doCut = 0;
	int preroll = 0, postroll = 0, doFloor = 0, ass_quality = 0, nify = 0; double damned = 0;
	if (!strcmp(argv[argc-2],"N-IFY"))
		nify = atoi(argv[argc-1]), argc -=2, 
		printf("N-ifying all bases of quality (%d).\n",nify);
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
	if (argc > 5)  // Well, the last thing they can do is define QUAL_RIGHT, right?
		tright = atoi(argv[5]);
	printf("The minimum length a read can be is %d bases.\n",minLen);
	if (!doCut) printf("Trimming while quality < %d on left, < %d on right.\n",tleft,tright);
	
	// Go through their files, C the contents. C what I did there? C?
	FILE *in = fopen(argv[1],"rb"), *out = fopen(argv[2],"wb");
	if (!in || !out) {puts("Can't open input or output file(s)!"); exit(1);}
	char *Head = malloc(UINT16_MAX), *Seq = malloc(UINT16_MAX), 
		*Shi7 = malloc(UINT16_MAX), *Qual = malloc(UINT16_MAX); 
	size_t crap = 0, decency = 0, totI = 0, totJ = 0, totLen = 0, length;
	int *Scores = calloc(doRoll,sizeof(*Scores));
	while (Head = fgets(Head,UINT16_MAX,in)) { 
		Seq = fgets(Seq, UINT16_MAX, in);
		Shi7 = fgets(Shi7, UINT16_MAX, in);
		Qual = fgets(Qual, UINT16_MAX, in);
		length = strlen(Qual);
		if (Qual[length-1]=='\n') --length;
		if (Qual[length-1]=='\r') --length;
		
		// Read--no, FEEL the quality (some of these methods have a lookahead).
		int i, j;
		if (doRoll) { // Do the rolling average thing. FILO queues, you know.
			int total = 0; int roll = 0, pos = 0;
			for (i = 0; i < length; i++) {
				Scores[pos] = Qual[i] - 33;
				total += Scores[pos];
				if (++pos == doRoll) pos = 0;
				if (roll == doRoll) total -= Scores[pos];
				else {++roll; if (!preroll) continue;}
				if ((float)total/roll > (float)tleft) break;
			}
			if (postroll) i -= roll - 1;
			if (i == length) {++crap; continue;}
			total = 0; roll = 0; pos = 0;
			for (j = length-1; j > i; j--) {
				Scores[pos] = Qual[j] - 33;
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
			for (i = 0; i < length; i++) {
				if (Qual[i] - 33 > tleft) ++row;
				else row = 0;
				if (row == doFloor) break;
			}
			i -= doFloor - 1;
			if (i == length) {++crap; continue;}
			row = 0;
			for (j = length-1; j > i; j--) {
				if (Qual[j] - 33 > tright) ++row;
				else row = 0;
				if (row == doFloor) break;
			}
			j += doFloor - 1;
		} 
		else if (doCut)  // Literally, idiotically slice off a specified number of bases.
			i = tright, j = length - tleft - 1;
		else { // By the way, the default trimming method is... go till ya hit somethin good.
			for (i = 0; i < length; i++)
				if (Qual[i] - 33 > tleft) break;
			if (i == length) {++crap; continue;}
			for (j = length-1; j > i; j--)
				if (Qual[j] - 33 > tright) break;
		}
		if (j == i || j - i < minLen) {++crap; continue;}
		++j;
		size_t abs = 0; for (int z = i; z < j; z++) abs += Qual[z];
		double thisAvgQ = (double)abs/(j - i);
		if (thisAvgQ - 33 < ass_quality) {++crap; continue;}
		damned += thisAvgQ;
		
		// Regurgitate the sequence. It's (prolly) OK enough if it made it this far down.
		fprintf(out,"%s",Head);
		for (int z = i; z < j; z++) fprintf(out,"%c",Qual[z] - 33 > nify ? Seq[z] : 'N');
		fprintf(out,"\n%s",Shi7);
		for (int z = i; z < j; z++) fprintf(out,"%c",Qual[z]);
		fprintf(out,"\n");
		++decency;
		totLen += j - i;
		totI += i, totJ += length - j;
	}
	printf("Done.\n%ld sequences were utter shi7, quality-wise.\n",crap);
	printf("And %ld sequences were decent.\n", decency);
	printf("The average trimmed length was %0.2f.\n", (double)totLen/decency);
	printf("On average, the cut bases were: %0.2f (left), %0.2f (right).\n", 
		(double)totI/decency, (double)totJ/decency);
	printf("The average read quality was %0.3f.\n-----\n", damned/decency - 33);
	free(Head); free(Seq); free(Shi7); free(Qual);
	return 0;
}
