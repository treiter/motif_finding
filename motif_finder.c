#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

static void generate_profile(int, int*, int**, char**, int);
static int profile_score(int**);


int k = 6;
int d = 0;
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
    printf("A: %d, C: %d, G: %d, T: %d, total: %d\n", a_count, c_count, g_count, t_count, (int)total_count);
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
    int sequenceToSkip = 0;
    generate_profile(sequenceToSkip, startIndices, profile, sequences, totalSequences);
    profile_score(profile);


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


static int profile_score(int** profile) {
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
	return -1;
}