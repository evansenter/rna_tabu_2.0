
#include "get_barrier.h"
#include "foldMod.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"

#define MAXLINE 3000

int disp_ptable2(short intable[500]){
  
  int x;
  short local[500];
  
  memcpy(local,intable,500*sizeof(short));
  
//  for (x=0;x<=local[0];x++)
  //  printf(" %d ",local[x]);

  for (x=1;x<=local[0];x++){
    if (local[x]==0)
      printf(".");
    else if (local[x]>x){
      printf("(");
    //  intable[local[x]]=1;
    }else
      printf(")");
  }
  printf("\n");

}

int main(int argc, char *argv[]){
  char sequence[500], start_str[500], end_str[500];
  int energy, e1, e2, barrier, tmp;
  eReturn eret;
  //for solution storage
  short route[1000][500];
  short best_route[1000][500];
  int best_energy=1000000;
  int route_len;
  int best_route_len;
  int best_k;

  //for best solution storage
  short best_route_all[1000][500];
  int best_energy_all=1000000;
  int best_route_len_all;
  int best_k_all;


  int k,i,t,j,l,a;
  int n;
  int print_iterations=1; // print all routes
  int reps=10; // default number of iterations
  int st; //start structure energy
  int tg; //end structure energy
  int k1=10;  // default minumum weight
  int k2=70;  // default maximum weight
  int th=10000; // default maximum energy barrier
  int user_dangles=2; // default dangling end treatment
  FILE *fp;
  char lines[4][MAXLINE], fname[100], comment[100];
  int fromFAA=0;
/* OLD INPUR FORMAT
  if (argc==8){
    strcpy(sequence, argv[1]);
    strcpy(start_str, argv[2]);
    n = strlen(sequence);
    energy_of_struct(sequence, start_str);//May initialize some things.
    strcpy(end_str, argv[3]);
    reps = atoi(argv[4]);
    k1= atoi(argv[5]);
    k2 = atoi(argv[6]);
    th = (int)(atof(argv[7])*100);
  } else if (argc>2 && strcmp(argv[1], "-f")==0){
    fromFAA=1;
    strcpy(fname, argv[2]);
    fp = fopen(fname, "r");
    for (i=0;i<4;i++){
     if (fgets(lines[i],MAXLINE,fp)==NULL){
	printf("Cannot open fasta file.\n  Check format, need start and stop structure after sequence.\n", fname);
	exit(1);
      }
    }
    strcpy(comment, lines[0]);
    sscanf(lines[1], "%s", sequence);
    sscanf(lines[2], "%s", start_str);
    sscanf(lines[3], "%s", end_str);
    n = strlen(sequence);
    energy_of_struct(sequence, start_str);//May initialize some things.
    if(argc>3){
      reps = atoi(argv[3]);
      k1= atoi(argv[4]);
      k2 = atoi(argv[5]);
      th = (int)(atof(argv[6])*100);
    }
  } else {
    printf("Usage:\n");
    printf(" %s sequence start_structure end_structure iterations lb_init_weight ub_init_weight energy_bound\n", argv[0]);
    printf("    (use parenthesis around structures)\n");
    printf(" or \n");
    printf(" %s -f fasta_file_with_start_end_structure iterations lb_init_weight ub_init_weight energy_bound\n", argv[0]);
    exit(1);
  }*/


  /** NEW INPUT FORMAT **/
  for (i = 1; i<argc; i++){
    if (argv[i][0] == '-'){
      if (strcmp(argv[i], "-sq") == 0) {
        if(++i<argc){
          strcpy(sequence,argv[i]);
        }
        if ((strcmp(sequence,"") == 0) || (sequence[0] == '-')) {
          printf("\nError in sequence!\n");
          usage(argv[0]);
        }
      }
      else if (strcmp(argv[i], "-ss") == 0) {
        if(++i<argc){
          strcpy(start_str,argv[i]);
        }
        if ((strcmp(start_str,"") == 0) || (start_str[0] == '-')) {
          printf("\nnError in start structure!\n");
          usage(argv[0]);
        }
      }
      else if (strcmp(argv[i], "-es") == 0) {
        if(++i<argc){
          strcpy(end_str,argv[i]);
        }
        if ((strcmp(start_str,"")==0) || (end_str[0] == '-')) {
          printf("\nError in end structure!\n");
          usage(argv[0]);
        }
      }
      else if (strcmp(argv[i], "-f") == 0){
        if(++i<argc){
          strcpy(fname,argv[i]);
        }
        fp = fopen(fname, "r");
        for (j=0;j<4;j++){
          if (fgets(lines[j],MAXLINE,fp)==NULL){
	    printf("Cannot open fasta file.\n  Check format, need start and stop structure after sequence.\n", fname);
            exit(1);
          }          
        }
        strcpy(comment, lines[0]);
        sscanf(lines[1], "%s", sequence);
        sscanf(lines[2], "%s", start_str);
        sscanf(lines[3], "%s", end_str);
        close(fp);
      }
      else if (strcmp(argv[i], "-d") == 0){
        if(++i<argc){
          if (sscanf(argv[i], "%d", &user_dangles)==0)
            usage(argv[0]);
        }
      }
      else if (strcmp(argv[i], "-it") == 0){
        if(++i<argc){
          if (sscanf(argv[i], "%d", &reps)==0)
            usage(argv[0]);
        }
      }
      else if (strcmp(argv[i], "-eb") == 0){
        if(++i<argc){
          if (sscanf(argv[i], "%d", &th)==0)
            usage(argv[0]);
          th*=100;
        }
      }
      else if (strcmp(argv[i], "-wmin") == 0){
        if(++i<argc){
          if (sscanf(argv[i], "%d", &k1)==0)
            usage(argv[0]);
        }
      }
      else if (strcmp(argv[i], "-wmax") == 0){
        if(++i<argc){
          if (sscanf(argv[i], "%d", &k2)==0)
            usage(argv[0]);
        }
      }
      else if (strcmp(argv[i], "-pi") == 0){
        if(++i<argc){
          if (sscanf(argv[i], "%d", &print_iterations)==0)
            usage(argv[0]);
        }
      }
      else {
         usage(argv[0]);
      }
    }
  }
  if ((strcmp(start_str,"")==0) || (strcmp(end_str,"")==0) || (strcmp(sequence,"")==0)){
    usage(argv[0]);
  }
//  printf ("%s %s %s %d\n",sequence, start_str, end_str, user_dangles);
  n = strlen(sequence);
//  energy_of_struct(sequence, start_str);//May initialize some things.

  /** END NEW INPUT FORMAT **/

  for(t=0;t<reps;t++){
     best_energy=1000000; 
     for(k=k1;k<k2;k+=2){  
      srand( (unsigned int) time( NULL )+t+k);
      //srand( (unsigned int) reps);
      //printf("Try for k = %d\n",k);  
      
      eret=getBarrierEnergy(sequence, start_str, end_str, route, &route_len,k,th,user_dangles);
      barrier=eret.max;
      st= eret.start;

      //printf(" Barrier1 = %d (max = %d start = %d\n", barrier,eret.max,eret.start);
      //for(i=0;i<route_len;i++)
      // disp_ptable2(route[i]);
      //getchar();
      if(barrier<best_energy){
         best_energy=barrier;
         best_route_len = route_len;
         best_k=k;
         for(i=0;i<route_len;i++)
           memcpy(best_route[i],route[i],(n+1)*sizeof(short));
      }
      eret=getBarrierEnergy(sequence, end_str, start_str, route, &route_len,k,th, user_dangles);
      barrier=eret.max;
      tg= eret.start;
      //printf(" Barrier2 = %d (max = %d start = %d\n", barrier,eret.max,st);
      //printf(" Barrier2 = %d \n", barrier);
      //for(i=0;i<route_len;i++)
      // disp_ptable2(route[i]);
      //getchar();
      if(barrier<best_energy){
         best_energy=barrier;
         best_route_len = route_len;
         best_k = k;
         for(i=0;i<route_len;i++)
           memcpy(best_route[i],route[route_len-i-1],(n+1)*sizeof(short));
      }
    }


    if(best_energy>50000)
      printf("No path found\n");
    else{
      if(print_iterations==1){
        if (fromFAA){
          printf("%s", comment);
        }
        printf("%s\n",sequence);
        for(i=0;i<best_route_len;i++)
          disp_ptable2(best_route[i]);
        printf("barrier is %5.2f -- route length = %d -- best k = %d th = %d -- real barrier is %5.2f\n", (float) best_energy/100.,best_route_len,best_k,th, (float) ((best_energy-st)/100.));
        //getchar();
      }

    }
    /* STORE BEST PATH OF ALL */
    if(best_energy<best_energy_all || (best_energy==best_energy_all && best_route_len <best_route_len_all)){
      best_energy_all = best_energy;
      for(i=0;i<best_route_len;i++)
        memcpy(best_route_all[i],best_route[best_route_len-i-1],(n+1)*sizeof(short));
      best_route_len_all=best_route_len;
      best_k_all=best_k;
    }
  }
  /* PRINT BEST PATH OF ALL */
  if(best_energy_all<=50000){
      printf("**BEST PATH**\n");
      printf("Start: %f\n",st/100.);
      printf("End: %f\n",tg/100.);
      if (fromFAA){
	printf("%s", comment);
      }
      printf("%s\n",sequence);
      for(i=0;i<best_route_len_all;i++)
       disp_ptable2(best_route_all[i]);
      printf("barrier is %5.2f -- route length = %d -- best k = %d th = %d -- real barrier is %5.2f\n", (float) best_energy_all/100.,best_route_len_all,best_k_all,th, (float) ((best_energy_all-st)/100.));
  }
  
  return 0;
}

void usage(char* name)
{

    printf("Usage:\n");
    printf(" %s -sq sequence -ss start_structure -es end_structure [-it iterations] [-wmin lb_init_weight] [-wmax ub_init_weight] [-eb energy_bound] [-d dangles] [-pi print_iterations]\n", name);
    printf("    (use quotes around structures)\n");
    printf(" or \n");
    printf(" %s -f fasta_file_with_comment_sequence_start_end_structures [-it iterations] [-wmin lb_init_weight] [-wmax ub_init_weight] [-eb energy_bound] [-d dangles] [-pi print_iterations]\n", name);
    printf(" default values for optional arguments: \n");
    printf("  iterations: 10 \n");
    printf("  lb_init_weight: 10 \n");
    printf("  ub_init_weight: 70 \n");
    printf("  energy_bound: 100 \n");
    printf("  print_iterations: 1 \n");
    printf("  dangles: 2 \n");
    exit(1);
}

