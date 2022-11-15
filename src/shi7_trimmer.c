#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#define SHOW_USAGE() {\
	printf( "\nSHort-read Iterative Trimmer (SHI7) v0.92d, by Gabe. Usage:\n");\
	printf( "shi7_trimmer in_seqs.fastq out_prefix MINLEN QUAL <QUAL_RIGHT> ...\n");\
	printf( "Choose one of the following optional trimming methods:\n");\
	printf( "   [ROLLING X [PREROLL/POSTROLL]] For rolling average over X bases\n" );\
	printf( "      PREROLL allows initial windows < X. POSTROLL adds back X to both sides.\n");\
	printf( "   [FLOOR X] Cut until X consecutive bases have minimum quality\n");\
	printf( "   [CUT] Treat 'QUAL' and 'QUAL_RIGHT' as number of bases to blindly cut\n");\
	printf( "[ASS_QUALITY X] Chuck reads under minimum average quality X (after trim).\n");\
	printf( "[N-IFY X] Magically transform crap quality (<= X) into ambiguous bases.\n");\
	printf( "[CASTN X] Magically transform Ns into quality score X.\n");\
	printf( "[DITCH X] Ditch all reads with ANY quality score below X.\n");\
	printf( "[STRIP] Strip the header names and replace with numerical indices.\n");\
	printf( "[HEAD X] Only output the first X good sequences.\n");\
	printf( "[AUTOADAP] Shear after common adapter sequences and tails [SE/PE].\n");\
	printf( "[ADAP ...] Shear after chosen adapter seq. Can be used up to 5x\n");\
	printf( "[ADAP2 ...] Shear after paired partial adapter cap. Can be used up to 5x\n");\
	printf( "[PALINCUT X] Shear after auto-detected paired-end adapter (Xtot-2end)\n");\
	printf( "[R2 in_seqsR2.fastq]: Invoke paired-end mode with specified R2.\n");\
	printf( "[OUTFASTA]: Output to FASTA format instead of FASTQ.\n");\
	printf("\nOh, you'll want to add 31 to all quality args if you're on ol' PHRED64.\n");\
	exit(1); \
}
char BUFFER[4096] = {0};
int main(int argc, char *argv[]) {
	// It's wise to start all your encounters with arguments.
	if (argc < 5 || argc > 25) SHOW_USAGE()
	int minLen = atoi(argv[3]), tleft = atoi(argv[4]), tright = tleft, doRoll = 0, doCut = 0;
	int preroll = 0, postroll = 0, doFloor = 0, nify = 0; double damned = 0., ass_quality = 0.;
	int head = INT32_MAX, strip = 0, ditch = 0, palincut = 0, pend=2, fasta_out = 0, doAdap2 = 0;
	char *pair = 0, *ZAdaps[5] = {0,0,0,0,0}, **Adaps = ZAdaps, 
		*AutoAdaps[5] = {"AGATCGGAAGAGC","CTGTCTCTTATA","TGGAATTCTCGG","AAAAAAAAAA","GGGGGGGGGG"};
	int AdapLens[5] = {13,12,12,10,10}, nAdaps2 = 0, doCastn = 0, castn = 0;
	if (!strcmp(argv[argc-1],"OUTFASTA"))
		fasta_out = 1, argc -=1, 
		printf("Outputting to FASTA files.\n");
	if (!strcmp(argv[argc-2],"R2"))
		pair = argv[argc-1], argc -= 2,
		printf("Pairing '%s' and '%s'.\n",argv[1],pair);
	if (!strcmp(argv[argc-2],"PALINCUT"))
		palincut = atoi(argv[argc-1]), argc -=2, 
		pend = palincut < 0 ? 1 : 2, palincut = abs(palincut), 
		printf("Auto-removing PE adapters of len >=%d or >=%d at end.\n",
			palincut,pend);
	if (!strcmp(argv[argc-1],"AUTOADAP"))
		Adaps = AutoAdaps, nAdaps2 = 5, argc -=1, 
		printf("Auto-removing some suspected adapters [SE/PE]. Ignoring further ADAP.\n");
	if (Adaps != AutoAdaps && nAdaps2 == 0) for (int i = 0; i < 5; ++i) 
		if (!strcmp(argv[argc-2],"ADAP") || !strcmp(argv[argc-2],"ADAP2"))
			doAdap2 |= !strcmp(argv[argc-2],"ADAP2"),
			Adaps[i] = argv[argc-1], AdapLens[i] = strlen(argv[argc-1]), ++nAdaps2, argc -=2,
			printf("Clipping %s adapter # %d seq '%s' of len %d.\n",
				doAdap2? "capped PE":"independent",i+1,Adaps[i],AdapLens[i]);
	if (!strcmp(argv[argc-2],"HEAD"))
		head = atoi(argv[argc-1]), argc -=2, 
		printf("Outputting only the first %d sequences.\n",head);
	if (!strcmp(argv[argc-1],"STRIP"))
		strip = 1, argc -=1, 
		printf("Stripping superflous header and + string contents.\n");
	if (!strcmp(argv[argc-2],"DITCH"))
		ditch = atoi(argv[argc-1]), argc -=2, 
		printf("Ditching ALL reads with ANY quality < %d.\n",ditch);
	if (!strcmp(argv[argc-2],"CASTN"))
		doCastn = 1, castn = atoi(argv[argc-1]); castn = castn < 0 ? 0 : castn, argc -=2, 
		printf("Casting N's to quality score %d.\n",castn);
	if (!strcmp(argv[argc-2],"N-IFY"))
		nify = atoi(argv[argc-1]), argc -=2, 
		printf("N-ifying all bases of quality <= %d.\n",nify);
	if (!strcmp(argv[argc-2],"ASS_QUALITY"))
		ass_quality = atof(argv[argc-1]), argc -=2, 
		printf("Chucking reads of ass quality (%g).\n",ass_quality);
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
	if (!pair) printf("Considering SE sample %s\n",argv[1]);
	printf("The minimum length a read can be is %d bases.\n",minLen);
	if (!doCut) printf("Trimming while quality < %d on left, < %d on right.\n",tleft,tright);
	
	// Go through their files, C the contents. C what I did there? C?
	sprintf(BUFFER,"%s%s.f%c",argv[2],pair? ".R1" : "",fasta_out? 'a':'q');
	gzFile in = gzopen(argv[1],"rb"), in2 = 0;
	FILE *out = fopen(BUFFER,"wb"), *out2 = 0;
	setvbuf(out,0,_IOFBF,1<<20);
	if (pair) 
		sprintf(BUFFER,"%s.R2.f%c",argv[2],fasta_out? 'a':'q'),
		in2=gzopen(pair,"rb"), out2 = fopen(BUFFER,"wb"),
		setvbuf(out2,0,_IOFBF,1<<20);
	if (!in || !out || (pair && (!in2 || !out2))) {
		puts("Can't open input or output file(s)!"); exit(1);}
	char *Head = malloc(UINT16_MAX), *Seq = malloc(UINT16_MAX), 
		*Shi7 = malloc(UINT16_MAX), *Qual = malloc(UINT16_MAX),
		*Head2 = 0, *Seq2 = 0, *Shi72 = 0, *Qual2 = 0; 
	if (pair) Head2 = malloc(UINT16_MAX), Seq2 = malloc(UINT16_MAX),
		Shi72 = malloc(UINT16_MAX), Qual2 = malloc(UINT16_MAX);
	if (palincut && !pair) puts("WARNING: Can't do palincut on SE reads!");
	long crap = 0, decency = 0, totI = 0, totJ = 0, totLen = 0, length=0, length2=LONG_MAX;
	int *Scores = calloc(doRoll,sizeof(*Scores));
	while (Head = gzgets(in,Head,UINT16_MAX)) { 
		Seq = gzgets(in,Seq, UINT16_MAX);
		Shi7 = gzgets(in,Shi7, UINT16_MAX);
		Qual = gzgets(in,Qual, UINT16_MAX);
		length = strlen(Qual);
		//printf("Seq: %s qual: %s len %d\n",Seq,Qual,length);
		if (Qual[length-1]=='\n') --length;
		if (Qual[length-1]=='\r') --length;
		//printf("--> Length is now %d\n",length);
		if (pair) {
			Head2 = gzgets(in2,Head2, UINT16_MAX);
			Seq2 = gzgets(in2,Seq2, UINT16_MAX);
			Shi72 = gzgets(in2,Shi72, UINT16_MAX);
			Qual2 = gzgets(in2,Qual2, UINT16_MAX);
			length2 = strlen(Qual2);
			if (Qual2[length2-1]=='\n') --length2;
			if (Qual2[length2-1]=='\r') --length2;
		}
		long L, R, L2=0, R2=0, adMin = INT32_MAX; 
		double darned = 0.; // A fitting euphamism
		
		// Prepass for the adapters only -- hey, we're adaptable.
		for (int r = 0; r <= (pair != 0); ++r) {
			char *TQual = Qual, *TSeq = Seq; long *TLen = &length; 
			if (r) TQual = Qual2, TSeq = Seq2, TLen = &length2;
			for (int z = 0; z < nAdaps2; ++z) { // Truncate at adaps
				char *adap = strstr(TSeq,Adaps[z]);
				long endpt = adap ? (adap-TSeq + 1) : *TLen;
				if (pair && doAdap2) {
					adMin = adMin < endpt ? adMin : endpt; // truncate both reads where first adap found
					// Truncate at partial matches up to 2 bases near the end
					if (!r && !adap) // Only run on R1; second read is redundant due to equality
						for (int w = *TLen - (AdapLens[z] < *TLen ? AdapLens[z] : 0) + 1; w < *TLen-1; ++w) {
							if (!memcmp(Adaps[z],TSeq+w,*TLen - w) && !memcmp(Seq+w,Seq2+w,*TLen-w))
								{adMin = adMin < w ? adMin : w; break;}
						//printf("--> adMin post-adaptrim? %ld\n",adMin);
					}
				}
				*TLen = *TLen < endpt ? *TLen : endpt;
			}
		}
		if (pair && palincut) { // prepass for PE palincut adapters
			long same = 0, minLen = length < length2 ? length : length2;
			minLen = minLen < adMin ? minLen : adMin;
			for (int i = 0; i < minLen; ++i) {
				same = (Seq[i]==Seq2[i]) ? (same + 1) : 0;
				if (same == palincut) {adMin = i-palincut+1; same = 0; break;}
			}
			if (same >= pend) adMin -= same;
		}
		length = adMin < length ? adMin : length;
		length2 = adMin < length2 ? adMin : length2;
		//printf("--> final pre-qc lens: %ld, %ld\n",length,length2);

		// Double the mates, half the fun? Quite the mechanism, a-pair-end-ly!
		for (int z = 0; z <= (pair != 0); ++z) {
			char *TQual = Qual, *TSeq = Seq; long TLen = length; 
			if (z) TQual = Qual2, TSeq = Seq2, TLen = length2;
			//TSeq[TLen] = 0, TQual[TLen] = 0;

			int nlPos = strchr(TSeq,'\n')-TSeq; // Truncate at embeded NLs
			TLen = TLen < nlPos ? TLen : nlPos;

			// The rains of cast-a-mere Ns
			if (doCastn) for (int i = 0; i < TLen; ++i) 
				TQual[i] = TSeq[i] == 'N' ? 33 + castn : TQual[i];

			// Read -- no, FEEL the quality (some have gotta 'lookahead').
			int i, j;
			if (doRoll) { // Do the rolling average thing. Fill-all FILO queues
				int total = 0; int roll = 0, pos = 0;
				int nlPos = strchr(Seq,'\n')-Seq;
				R = R < nlPos ? R : nlPos;
				for (i = 0; i < TLen; i++) {
					Scores[pos] = TQual[i] - 33;
					total += Scores[pos];
					if (++pos == doRoll) pos = 0;
					if (roll == doRoll) total -= Scores[pos];
					else {++roll; if (!preroll) continue;}
					if ((float)total/roll > (float)tleft) break;
				}
				if (postroll) i -= roll;
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
				if (postroll) j += roll;
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
				//printf("--> post-floor i: %ld, j: %ld\n",i,j);
			} 
			else if (doCut)  // Idiotically slice off a set number of bases.
				i = tright, j = TLen - tleft - 1;
			else { // The default trimming method is 'go till ya hit somethin good'
				for (i = 0; i < TLen; i++)
					if (TQual[i] - 33 > tleft) break;
				if (i == TLen) {++crap; goto HELL;}
				for (j = TLen-1; j > i; j--)
					if (TQual[j] - 33 > tright) break;
			}
			if (j == i || j - i < minLen || j + 1 > TLen) {++crap; goto HELL;}
			++j;
			long abs = 0; for (int z = i; z < j; z++) abs += TQual[z];
			double thisAvgQ = (double)abs/(j - i);
			if (thisAvgQ - 33 < ass_quality) {++crap; goto HELL;}
			if (ditch) for (int z = i; z < j; ++z) 
				if (TQual[z] - 33 < ditch) {++crap; goto HELL; }
			darned += thisAvgQ;
			if (z) L2 = i, R2 = j; else L = i, R = j;
		}
		damned += darned; 
		
		// Regurgitate the sequence(s). They're (prolly) OK enough if they made it this far down.
		for (int z = L; z < R; z++) Seq[z] = Qual[z] - 33 > nify ? Seq[z] : 'N';
		Seq[R]=0; Qual[R]='\n', Qual[R+1] = 0; 
		//printf(" ------> seq [%d] R1 '%s', L = %d, R = %d, len %d\n",crap+decency,Seq,L,R, R-L);
		if (fasta_out) Head[L]='>', Shi7[0]=0, Qual[L]=0;
		if (strip) fprintf(out,"%c%d\n%s\n%s%s",fasta_out? '>':'@',crap+decency,
			Seq + L,fasta_out? "":"+\n",Qual + L);
		else fprintf(out,"%s%s\n%s%s",Head,Seq + L,Shi7,Qual + L);
		if (pair) {
			for (int z = L2; z < R2; z++) Seq2[z] = Qual2[z] - 33 > nify ? Seq2[z] : 'N';
			Seq2[R2]=0; Qual2[R2]='\n'; Qual2[R2+1]=0;
			if (fasta_out) Head2[L2]='>', Shi72[0]=0, Qual2[L2]=0;
			if (strip) fprintf(out2,"%c%d\n%s\n%s%s",fasta_out? '>':'@',crap+decency,
				Seq2 + L2,fasta_out? "":"+\n",Qual2 + L2);
			else fprintf(out2,"%s%s\n%s%s",Head2,Seq2 + L2,Shi72,Qual2 + L2);
		}
		++decency;
		totLen += (R - L) + (R2 - L2);
		totI += L+L2, totJ += (length - R) + (length2 - R2);

		HELL:NULL; // HEREIN LIE DEVILS
		if (decency >= head) break;
	}
	printf("Done. %ld sequences were utter shi7, quality-wise.\n",crap);
	printf("And %ld sequences were decent.\n", decency);
	printf("The average trimmed length was %0.2f.\n", (double)totLen/decency/(1+(pair!=0)));
	printf("On average, the cut bases were: %0.2f (left), %0.2f (right).\n", 
		(double)totI/decency/(1+(pair!=0)), (double)totJ/decency/(1+(pair!=0)));
	printf("The average read quality was %0.3f.\n-----\n", damned/decency/(1+(pair!=0)) - 33);
	return 0;
}