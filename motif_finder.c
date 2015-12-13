#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <sys/time.h>
#include <stdbool.h>

//helper functions
static void generate_profile(int*, int**, char**, int);
static void remove_from_profile(int, int*, int**, char**);
static void add_to_profile(int, int*, int**, char**);
static double profile_score(int**, int, double*, int*);
static bool isValInArray(int, int*);

//default global values
int k = 6;
int d = 0;
int t = 2;

int main(int argc, char** argv) {
    //structs used for timing
    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    FILE * fp;

    int option = 0;
    //parse arguments
    while((option = getopt(argc, argv, "k:d:t:")) != -1) {
        switch(option) {
            case 'k': 
                k = atoi(optarg);
                break;
            case 'd':
                d = atoi(optarg);
                break;
            case 't':
                t = atoi(optarg);
                break;
            default:
                printf("Usage: motif_finder [-k motif_size] [-d wild_cards] [-t max_time] input_file");
                break;
        }
    }

    if((fp = fopen(argv[argc-1], "r")) == NULL) {
        printf("FILE: \"%s\" DOESN'T EXIST\n", argv[argc - 1]);
        exit(1); //Failed to open file.
    }

    //count how big our buffer needs to be (every other line)
    char line[80];
    int totalSequences = 0;
    while(fgets(line, 80, fp) != NULL) {
    	if(line[0] == '>') {
    		totalSequences++;
    	}
    }
    rewind(fp);//reset fp to start of file
    //read in the sequences. account for multi-line sequences
    char* sequences[totalSequences];
    int sequenceNum = 0;
    int sequenceLength = 0;
    int beginningOfSequence = 1;
    while(fgets(line, 80, fp) != NULL) {
    	//sequence header
    	if(line[0] == '>') {
    		sscanf(line, ">Sequence%d length %d", &sequenceNum, &sequenceLength);
    		sequenceNum--; //keep 0 indexed
            //note: sequenceLength is multiplied by 2, because it seemed that some of the "lengths" were off in the file
    		char* temp  = malloc((sequenceLength*2)*sizeof(char));
            if(temp == NULL) {
                return -1;
            }
            sequences[sequenceNum] = temp;
    		beginningOfSequence = 1;
        //read the line (if its the start of a sequence)
    	} else if(beginningOfSequence) {
    		strcpy(sequences[sequenceNum], line);
    		beginningOfSequence = 0;
        //otherwise append it to the current sequences
    	} else {
    		strcat(sequences[sequenceNum], line);
    	}
        //clear line (to be safe)
        memset(line, 0, 80*sizeof(char));
    }
    fclose(fp);

    //determine frequency
    int a_count = 0, c_count = 0, g_count = 0, t_count = 0;
    double total_count = 0; //to force double division later for probabilities
    for(int i = 0; i < totalSequences; i++) {
    	for(int j = 0; j < strlen(sequences[i]); j++) {
    		switch(sequences[i][j]) {
    			case 'A':
    				a_count++;
    				total_count++;
    				break;
    			case 'C':
    				c_count++;
    				total_count++;
    				break;
    			case 'G':
    				g_count++;
    				total_count++;
    				break;
    			case 'T':
    				t_count++;
    				total_count++;
    				break;
    		}
    	}
    }
    //only calculate this once
    double a_freq = a_count/total_count;
    double c_freq = c_count/total_count;
    double g_freq = g_count/total_count;
    double t_freq = t_count/total_count;
    double freqs[4];
    freqs[0] = a_freq;
    freqs[1] = c_freq;
    freqs[2] = g_freq;
    freqs[3] = t_freq;

    //initialize our profile matrix
    int* profile[4];
    for(int i = 0; i < 4; i++) {
    	profile[i] = (int*)malloc((k+1)*sizeof(int));//slot 0 will be the total count
    	memset(profile[i], 0, (k+1)*sizeof(int));
    }

    //pick random starting locations
    int startIndices[totalSequences];
    srand(time(NULL));
    for(int i = 0; i < totalSequences; i++) {
    	startIndices[i] = rand() % (strlen(sequences[i])-k);
    }
    generate_profile(startIndices, profile, sequences, totalSequences);
    //kind of hacky to avoid errors when d == 0
    int* dontCares = malloc(((d==0)?1:d)*sizeof(int));
    int* bestDontCares = malloc(((d==0)?1:d)*sizeof(int));
    //grab the inital score
    double logScore = profile_score(profile, totalSequences, freqs, dontCares);
    double bestLogScore = logScore;
    //more timing
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);

    //proceed until we're out of time
    while(tval_result.tv_sec < t) {
        //for each sequence
        for(int i = 0; i < totalSequences; i++) {
            int oldStartIndex = startIndices[i];
            int currBestStart = oldStartIndex;
            //try each starting position
            remove_from_profile(i, startIndices, profile, sequences);
            for(int j = 0; j < strlen(sequences[i])-k; j++) {
                if(j != oldStartIndex) {
                    startIndices[i] = j;
                    add_to_profile(i, startIndices, profile, sequences);
                    logScore = profile_score(profile, totalSequences, freqs, dontCares);
                    //keep track of the best score so far
                    if(logScore > bestLogScore) {
                        currBestStart = j;
                        bestLogScore = logScore;
                        memcpy(bestDontCares, dontCares, d*sizeof(int));
                    }
                    remove_from_profile(i, startIndices, profile, sequences);
                }
            }
            //reset profile with the best new starting index for sequences i
            startIndices[i] = currBestStart;
            add_to_profile(i, startIndices, profile, sequences);
            //check if we're out of time yet
            gettimeofday(&tval_after, NULL);
            timersub(&tval_after, &tval_before, &tval_result);
            if(tval_result.tv_sec >= t) {
                break;
            }
        }
    }
    //done with algorithm, now just print stuff
    //figure out the motif
    char* motif = (char*)malloc((k+1)*sizeof(char));
    if(motif == NULL) {
        return -1;
    }
    for(int i = 0; i < k; i++) {
        int biggestChar = 0;
        motif[i] = 'A';
        if(isValInArray(i, bestDontCares)) {
            motif[i] = '*';
            continue;
        }
        if(profile[1][i+1] > profile[biggestChar][i+1]) {
            biggestChar = 1;
            motif[i] = 'C';
        }
        if(profile[2][i+1] > profile[biggestChar][i+1]) {
            biggestChar = 2;
            motif[i] = 'G';
        }
        if(profile[3][i+1] > profile[biggestChar][i+1]) {
            biggestChar = 3;
            motif[i] = 'T';
        }
    }
    //print results
    motif[k] = '\0';
    printf("Best motif of length %d with %d don't cares is %s\n", k, d, motif);
    printf("Log likelihood is %lf\n", bestLogScore);
    printf("Loci of the best motif are here:\n");

    for(int i = 0; i < totalSequences; i++) {
        printf("%d\n", startIndices[i]);
        free(sequences[i]);
    }
    //clean up heap
    for(int i = 0; i < 4; i++) {
        free(profile[i]);
    }
    free(dontCares);
    free(bestDontCares);
    free(motif);
    return 0;
}

static bool isValInArray(int val, int* array) {
    for(int i = 0; i < d; i++) {
        if(array[i] == val) {
            return true;
        }
    }
    return false;
}

static void generate_profile(int* startIndices, int** profile, char** allSequences, int numSequences) {
	for(int i = 0; i < numSequences; i++) {
		for(int j = 0; j < k; j++) {
			char currChar = allSequences[i][j+startIndices[i]];
			switch(currChar) {
				case 'A':
    				profile[0][j+1]++;
    				profile[0][0]++;
    				break;
    			case 'C':
    				profile[1][j+1]++;
    				profile[1][0]++;
    				break;
    			case 'G':
    				profile[2][j+1]++;
    				profile[2][0]++;
    				break;
    			case 'T':
    				profile[3][j+1]++;
    				profile[3][0]++;
    				break;
			}
		}
		
	}
}

//This should only be called if the profile doesn't include information from the sequenceToAdd
static void add_to_profile(int sequenceToAdd, int* startIndices, int** profile, char** allSequences) {
	for(int i = 0; i < k; i++) {
		char currChar = allSequences[sequenceToAdd][i+startIndices[sequenceToAdd]];
		switch(currChar) {
			case 'A':
				profile[0][i+1]++;
				profile[0][0]++;
				break;
			case 'C':
				profile[1][i+1]++;
				profile[1][0]++;
				break;
			case 'G':
				profile[2][i+1]++;
				profile[2][0]++;
				break;
			case 'T':
				profile[3][i+1]++;
				profile[3][0]++;
				break;
		}
	}
}

//This should only be called if the profile already includes information from the sequenceToRemove
static void remove_from_profile(int sequenceToRemove, int* startIndices, int** profile, char** allSequences) {
	for(int i = 0; i < k; i++) {
		char currChar = allSequences[sequenceToRemove][i+startIndices[sequenceToRemove]];
		switch(currChar) {
			case 'A':
				profile[0][i+1]--;
				profile[0][0]--;
				break;
			case 'C':
				profile[1][i+1]--;
				profile[1][0]--;
				break;
			case 'G':
				profile[2][i+1]--;
				profile[2][0]--;
				break;
			case 'T':
				profile[3][i+1]--;
				profile[3][0]--;
				break;
		}
	}
}


static double profile_score(int** profile, int numSequences, double* freqs, int* dontCares) {
    double totalScore = 1;
    double scoreArray[k];
    double numSequencesAsDouble = (double)numSequences;
    for(int i = 0; i < k; i++) {
        double bestTotalProb = profile[0][i+1]/numSequencesAsDouble;
        bestTotalProb /= freqs[0];
        double prob1 = profile[1][i+1]/numSequencesAsDouble/freqs[1];
        double prob2 = profile[2][i+2]/numSequencesAsDouble/freqs[2];
        double prob3 = profile[3][i+3]/numSequencesAsDouble/freqs[3];
        if(prob1 > bestTotalProb) {
            bestTotalProb = prob1;
        }
        if(prob2 > bestTotalProb) {
            bestTotalProb = prob2;
        }
        if(prob3 > bestTotalProb) {
            bestTotalProb = prob3;
        }
        totalScore *= bestTotalProb;
        scoreArray[i] = bestTotalProb;
    }
    //figure out don't cares
    for(int i = 0; i < d; i++) {
        //get the minimum score in the array
        double minScore = INT_MAX;
        //dont include the edges
        for(int j = 1; j < k-1; j++) {
            if(scoreArray[j] < minScore) {
                minScore = scoreArray[j];
                dontCares[i] = j;
            }
        }
        //remove the don't cares from the score
        totalScore /= scoreArray[dontCares[i]];
        //don't consider this in the future
        scoreArray[dontCares[i]] = INT_MAX;
    }
	return log2(totalScore);
}