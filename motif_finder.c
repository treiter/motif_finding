#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

int main(int argc, char** argv) {

    //Set up defaults
    int k = 6;
    int d = 0;
    int t = 2;

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
        printf("FILE: %s DOESN'T EXIST\n", argv[argc - 1]);
    }


    fclose(fp);
    return 0;
}