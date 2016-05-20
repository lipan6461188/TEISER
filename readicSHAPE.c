#include "stdio.h"
#include "stdlib.h"

#include "structures.h"
#include "nucleotides.h"
#include "dataio.h"
#include "sequences.h"
#include "statistics.h"
#include "readicSHAPE.h"


int search_and_import_Seq(s_sequence **sequences, int t_seq_count, char *geneName, float *icValues, int icGeneLen, s_sequence *tmp);

/* t_seq_count is the number of sequences */
int read_icSHAPE( char* icshapeFile, s_sequence **sequences, int t_seq_count) {
    char* buff ;
    FILE  *fp ;
    int mynmax = 100000;
    char* genename ;
    float* icSHAPE_Values;
    int gene_length;
    int valid_ic_nums=0;   /* the amount of icshape which there is a fasta seq match it */
    
    /*Initiaze The Buff*/
    buff = (char*)malloc(mynmax * sizeof(char)) ;
    
    printf("icSHAPE file: %s\n",icshapeFile);
    fp = fopen(icshapeFile, "r") ;
    if (!fp){
        printf("read_icSHAPE: please enter a valid filename (%s invalid)\n", icshapeFile) ;
        exit(0) ;
    }
    
    s_sequence *s=NULL;
    
    for (int i = 0; i < 5062; i++) {
        //printf("\n%s\n",sequences[i]->name);
        if ( !strcmp(sequences[i]->name, "ENSMUST00000092033"))
        {
            s = sequences[i];
            break;
        }
    }
    
    
    while (!feof(fp)){
        //Read a new line
        fgets(buff, mynmax, fp);
        if (feof(fp))
            break;

        chomp(buff);
        genename = mystrtok(buff, '\t');
        char *tok = mystrtok(NULL, '\t');
        gene_length = 0;
        
        /*Initiaze The icSHAPE Reactivity Values*/
        icSHAPE_Values = (float*)malloc((int)(mynmax/sizeof(float)) * sizeof(float));
        
        while (tok != NULL)
        {
         //   printf("read in a new value: %s\n",tok);
            /* return False when equal*/
            if (strcmp(tok, "NULL")) {
                //tok is not NULL
                icSHAPE_Values[gene_length] = atof(tok);
            }
            else{
                //tok is NULL
                icSHAPE_Values[gene_length] = -1;
            }
            gene_length++;
           // printf("free tock\n");
            free(tok);
            tok = mystrtok(NULL, '\t');
        }
        if( search_and_import_Seq(sequences,t_seq_count,genename,icSHAPE_Values,gene_length,s) )
        {
            valid_ic_nums++;
 /*
            if ( !strcmp(genename, "ENSMUST00000092033"))
            {
                printf("\n%s:\t",genename);
                for (int j=0; j<gene_length;j++)
                {
                    printf("%.2f \t",s->icSHAPE[j]);
                }
                printf("\n");
                
                if( s->icSHAPE == icSHAPE_Values)
                {
                    printf("\n EQUAL \n");
                }
            }
            printf("%s %.2f \n",genename,s->icSHAPE[0]);
*/
        }
        else
            free(icSHAPE_Values);
        //printf("free geneName\n");
        free(genename);
    }
    
  /*

            printf("third time\n%s:\t",s->name);
            for (int j=0; j<s->length;j++)
            {
                printf("%.2f \t",s->icSHAPE[j]);
            }
            printf("\n");
    
    */
    
    fclose(fp);
    return valid_ic_nums;
}

int search_and_import_Seq(s_sequence **sequences, int t_seq_count, char *geneName, float *icValues, int icGeneLen, s_sequence *tmp)
{
    char *seqGeneName;
    int i;
    /*First, find the sequence which match the icSHAPE*/
    for (i=0 ; i< t_seq_count ; i++){
        seqGeneName = strdup(sequences[i]->name) ;
        if (!strcmp(seqGeneName, geneName)){
            
            
            if( tmp == sequences[i] )
            {
                printf("Find: %s\n", sequences[i]->name);
            }
            
         //   printf("free seqGeneName\n");
            free(seqGeneName);
                /* Find It!*/
            /* If the length is not equal then exit */
            if( sequences[i]->length != icGeneLen )
            {
                printf("read_icSHAPE: The Length of Gene %s: length not equal:%d & %d\n",geneName,sequences[i]->length,icGeneLen);
                exit(0);
            }
            /*Now, Import the icSHAPE into Sequence*/
            free(sequences[i]->icSHAPE);
           // printf("load: %s\n", geneName);
            sequences[i]->icSHAPE = icValues;
            /*
            if ( !strcmp(geneName, "ENSMUST00000092033"))
            {
                printf("\n%s: \t",geneName);
                for (int j=0; j<icGeneLen;j++)
                {
                    printf("%.2f \t",icValues[j]);
                }
                printf("\n");
            }
            */
            break;
        }
    }
    if(i == t_seq_count)
        return 0;
    else
        return 1;
    
}














