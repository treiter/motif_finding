#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <math.h>

static void generate_profile(int, int*, int**, char**, int);
static void remove_from_profile(int, int*, int**, char**);
static void add_to_profile(int, int*, int**, char**);
static double profile_score(int**, int, double*, int*);


int k = 6;
int d = 2;
int t = 2;

int main(int argc, char** argv) {

    //Set up defaults


    FILE * fp;

    int option = 0;

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
    		sequences[sequenceNum] = malloc((sequenceLength+1)*sizeof(char));
    		beginningOfSequence = 1;
    	} else if(beginningOfSequence) {
    		strcpy(sequences[sequenceNum], line);
    		beginningOfSequence = 0;
    	} else {
    		strcat(sequences[sequenceNum], line);
    	}
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

    // printf("A: %d, C: %d, G: %d, T: %d, total: %d\n", a_count, c_count, g_count, t_count, (int)total_count);
    // printf("A: %lf, C: %lf, G: %lf, T: %lf\n", a_freq, c_freq, g_freq, t_freq);
    //initialize our profile matrix
    int* profile[4];
    for(int i = 0; i < 4; i++) {
    	profile[i] = (int*)malloc((k+1)*sizeof(int));//slot 0 will be the total count
    	memset(profile[i], 0, (k+1)*sizeof(int));
    }

    int startIndices[totalSequences];
    srand(time(NULL));
    for(int i = 0; i < totalSequences; i++) {
    	startIndices[i] = rand() % (strlen(sequences[i])-k);
    }
    //just testing from here out
    int sequenceToSkip = 0;
    generate_profile(sequenceToSkip, startIndices, profile, sequences, totalSequences);
    //following 2 lines are for testing
    remove_from_profile(1, startIndices, profile, sequences);
    add_to_profile(1, startIndices, profile, sequences);
    int dontCares[(d==0)?1:d];//kind of hacky to avoid errors
    int* bestDontCares = malloc(((d==0)?1:d)*sizeof(int));
    double bestLogScore = 0;
    double logScore = profile_score(profile, totalSequences, freqs, dontCares);
    printf("Score: %lf\n", logScore);
    printf("Dont cares:");
    for(int i = 0; i < d; i++) {
        printf(" %d", dontCares[i]);
    }
    printf("\n");

    //main loop goes something like while(end_condition)
    //  for(each sequence)
    //      align sequence with existing profile for best score
    //          ex: check each possible start index, compare score to bestScore
    //              keep track of "best" dontcares as well (so we dont lose that info)



    return 0;
}



static void generate_profile(int sequenceToSkip, int* startIndices, int** profile, char** allSequences, int numSequences) {
	for(int i = 0; i < numSequences; i++) {
		if(i != sequenceToSkip) {
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
        int biggestChar = 0;
        if(profile[1][i+1] > profile[biggestChar][i+1]) {
            biggestChar = 1;
        }
        if(profile[2][i+1] > profile[biggestChar][i+1]) {
            biggestChar = 2;
        }
        if(profile[3][i+1] > profile[biggestChar][i+1]) {
            biggestChar = 3;
        }
        double piProb = profile[biggestChar][i+1]/numSequencesAsDouble;
        double totalProb = piProb/freqs[biggestChar];//maybe want the max of this, not of charColumnCount
        totalScore *= totalProb;
        scoreArray[i] = totalProb;
        printf("%lf\n", totalProb);
    }
    //figure out don't cares
    for(int i = 0; i < d; i++) {
        //get the minimum score in the array
        double minScore = INT_MAX;
        //dont include the edges
        for(int j = 1; j < k-1; j++) {
            if(scoreArray[j] < minScore) {
                printf("%d, %d\n", i, j);
                minScore = scoreArray[j];
                dontCares[i] = j;
            }
        }
        //remove the don't cares from the score
        totalScore /= scoreArray[dontCares[i]];
        //don't consider this in the future
        scoreArray[dontCares[i]] = INT_MAX;
    }



	//JUST PRINTS OUT MATRIX
	for(int i = 0; i < 4; i++) {
		switch(i) {
			case 0:
				printf("A:\t");
				break;
			case 1:
				printf("C:\t");
				break;
			case 2:
				printf("G:\t");
				break;
			case 3:
				printf("T:\t");
				break;
		}
		for(int j = 0; j < k+1; j++) {
			printf("%d\t", profile[i][j]);
		}
		printf("\n");
	}


	//TODO: score matrix.
	return log2(totalScore);
}