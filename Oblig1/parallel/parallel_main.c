#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>


typedef struct {

  float** image_data;
  char* storage; // Brukes for kontinuerlig lagring
  int m; // height
  int n; // length
} image;


void import_JPEG_file(const char *filename, unsigned char **image_chars, int *image_height,
   int *image_width, int *num_components);

void export_JPEG_file(const char *filename, unsigned char *image_chars, int image_height,
   int image_width, int num_components, int quality);



void allocate_image(image *u, int m, int n);
void deallocate_image(image *u);
void convert_jpeg_to_image(const unsigned char* image_chars, image *u);
void convert_image_to_jpeg(const image *u, unsigned char* image_chars);
void iso_diffusion_denoising(image *u, image *u_bar, float kappa);


int main (int argc , char *argv[]) {

  int m, n, c, iters;
  int rest;
  int my_m, my_n, my_rank, num_procs;
  int *children_data;  //placeholder for information on size of childrens data sets
  int *extra_data; // plasssholder for informasjon om de ekstra linjene hver prosess skal ha
  int counter;
  int proc_m; // Proc 0s way to caluclate m for each process
  int extra ; // is the number of extra lines a prosess gets
  float kappa;
  image u, u_bar, whole_image;
  unsigned char *image_chars;
  unsigned char *my_image_chars;
  char *input_jpeg_filename;
  char *output_jpeg_filename;



  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  MPI_Status status;


  /*
    Checking if input paramethers are correct and of sufficent size
  */

  if(argc != 5){
    fprintf(stderr, "Usage: arguments %s\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  if((kappa = atof(argv[1])) <= 0){
    fprintf(stderr, "Usage: kappa too small or invalid : '%s'\n",argv[1]);
    exit(EXIT_FAILURE);
  }

  if((iters = atoi(argv[2])) < 3){
    fprintf(stderr, "Usage: iter too small must be run with at least 3 processes: '%s'\n",argv[2]);
    exit(EXIT_FAILURE);
  }

  input_jpeg_filename = argv[3];

  output_jpeg_filename = argv[4];


  /*
    All child process are ready to recieve their m size and ther extra size.
    Using this to send and return the right ammount of data between children and parent
  */

  if(my_rank > 0){
    MPI_Recv(&my_m, 1, MPI_INT, 0, 666, MPI_COMM_WORLD, &status);
    MPI_Recv(&extra, 1, MPI_INT, 0, 123, MPI_COMM_WORLD, &status);

  }


  if(my_rank == 0){  // Inne i process 0


    import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
    allocate_image(&whole_image, m, n);

/*
 * Plassholder for lengden på alle barns data.
 */
    if((children_data =  malloc(num_procs* sizeof(int))) == NULL){
      perror("Failed to allocate sizes chart for each process data size");
      exit(EXIT_FAILURE);
    }


    if((extra_data =  malloc(num_procs* sizeof(int))) == NULL){
      perror("Failed to allocate sizes chart for each process data size");
      exit(EXIT_FAILURE);
    }


    rest = m%(num_procs-1);  //-1 siden proc 0 ikke skal regne


  /*
   *  Dele opp data til barneprosessene og legge inn i arrayet for oversikt til senere
   */


    for(counter = 1;counter < num_procs; counter ++){

      if(counter <= rest ){ // hvis man er mindre eller lik  mudolo så gir jeg de en ekstra linje for rettferdighet

          proc_m = (m/(num_procs-1)) +1;

          if(counter == 1 || counter == num_procs-1){  // er man borderline skal men ha +1 til
            proc_m += 1;
            extra += 1;

          }

          else{ // e r man i midten skal man ha +2 for overlapping
            proc_m += 2;
            extra +=2;

          }

          MPI_Send(&proc_m, 1, MPI_INT, counter, 666, MPI_COMM_WORLD);
          children_data[counter] = proc_m;
          MPI_Send(&extra, 1, MPI_INT, counter, 123, MPI_COMM_WORLD);  // sender hvor mange de har fått ekstra for overlapping mellom prosessene
          extra_data[counter] = extra;
          extra = 0;

      } // end for alle som skal ha ekstra pga modulo og de er sjekket for å være borderline

      else{  // om man er storre en modulo skal man ha vanlig antall + de ma man evnt trenger for overlapp

          proc_m = (m/(num_procs-1));

          if(counter == 1 || counter == num_procs-1){  // er man borderline skal men ha +1
            proc_m += 1;
            extra += 1;

          }
          else{ // er man i midten skal man ha +2

            proc_m += 2;
            extra +=2;
          }

          MPI_Send(&proc_m, 1, MPI_INT, counter, 666, MPI_COMM_WORLD);   //sender storrelsen de skal allokere
          children_data[counter] = proc_m;
          MPI_Send(&extra, 1, MPI_INT, counter, 123, MPI_COMM_WORLD);  // sender hvor mange de har fått ekstra for overlapping mellom prosessene
          extra_data[counter] = extra;
          extra = 0; // resetting for next iteration if for-loop

      }

    } // end for-loop

  } // ute av process 0


  MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);


  if(my_rank > 0){ //1---x


    my_n = n;

    /*
    *  Allokering av plass for delbildet i arbeiderene
    */

    allocate_image (&u, my_m, my_n);
    allocate_image (&u_bar, my_m, my_n);


    /*
     * Allokering av plass i et kontonuerlig array
     */

     if((my_image_chars = (unsigned char*)malloc((my_m*my_n)*sizeof(unsigned char))) == NULL){
      perror("Failed malloc pointer to matrix");
      exit(EXIT_FAILURE);
    }

    /*
     * Recieve all image chars for iso_diffusion_denoising.
     */

     MPI_Recv(&my_image_chars[0], my_n*my_m, MPI_UNSIGNED_CHAR, 0,
              420, MPI_COMM_WORLD, &status);


      convert_jpeg_to_image(my_image_chars, &u);

       int upper = 1;
       int bottom = num_procs-1;

       u.m = my_m;
     	 u.n = my_n;

       for(counter = 0; counter < iters; counter++){

         iso_diffusion_denoising(&u, &u_bar, kappa);

         if(my_rank != upper){

           // Her sender vi oppover, noe vi ikke kan gjøre for den øveste prosessen

           MPI_Sendrecv(&u.image_data[1][0], my_n*4, MPI_UNSIGNED_CHAR, my_rank-1, 666,
             &u.image_data[0][0], my_n*4, MPI_UNSIGNED_CHAR, my_rank-1 , 666 , MPI_COMM_WORLD, &status);
         }

         if(my_rank != bottom){

           //Her sender vi nedover, noe som vi ikke kan gjøre i den nederste prosessen

            MPI_Sendrecv(&u.image_data[my_m-2][0], my_n*4,MPI_UNSIGNED_CHAR, my_rank+1, 666,
               &u.image_data[my_m-1][0],my_n*4, MPI_UNSIGNED_CHAR, my_rank+1 , 666 , MPI_COMM_WORLD, &status);
         }

       }//Ute av ISOdenoising


      /*
       *   Sending result back to process 0
       */

    	convert_image_to_jpeg(&u, my_image_chars);

      if(my_rank == upper){
        MPI_Send(&my_image_chars[0], ((my_m-extra)*my_n), MPI_UNSIGNED_CHAR,0,1337,MPI_COMM_WORLD);    //Kan starte på pos 0, og sende til
      }

      else{
        MPI_Send(&my_image_chars[my_n], ((my_m-extra)*my_n), MPI_UNSIGNED_CHAR,0,1337,MPI_COMM_WORLD);
      }


  }  // Ute av children


  else{ // In prosess 0



    counter = 0;
    int slave;

    /*
     *  Process 0 sends the correct ammount of data to the children
     *  counter += ((children_data[slave]-2)*n) is for setting the counter to the correct place
     *  so we are sure overlapping data is included in in the processes
     */

    for(slave = 1; slave <  num_procs; slave++){

      MPI_Send(&image_chars[counter], children_data[slave]*n, MPI_UNSIGNED_CHAR, slave,
                420, MPI_COMM_WORLD);

      counter += ((children_data[slave]-2)*n);

    }


    /*
     * get data back from children
     */
	counter = 0;

  	 for(slave = 1; slave <  num_procs; slave++){

        MPI_Recv(&image_chars[counter], (children_data[slave] - extra_data[slave])*n, MPI_UNSIGNED_CHAR, slave, 1337, MPI_COMM_WORLD, &status);
        counter += (((children_data[slave])-2)*n);

      }

    	export_JPEG_file(output_jpeg_filename, image_chars, m, n, 1, 75);

 /*
  * Deallocating image and freeing the arrays containing data
  */
      deallocate_image (&whole_image);
      free(children_data);
      free(extra_data);

  }  //Out of process 0


  /* Dealloc in child processes */
  if(my_rank > 0){
    free(my_image_chars);
    deallocate_image (&u);
    deallocate_image (&u_bar);
  }



  MPI_Finalize();

  return 0;
}


 void allocate_image(image *u, int m, int n){


   int i;
   memcpy(&u->m,&m,sizeof(int));
   memcpy(&u->n,&n,sizeof(int));


   if((u->image_data = (float**)malloc(m * sizeof(float*))) == NULL){
    perror("Failed malloc pointer to matrix");
    exit(EXIT_FAILURE);
   }


   for(i = 0; i < m; i++){
     if((u->image_data[i] = (float*)malloc(n * sizeof(float))) == NULL){
       perror("Inner malloc");
       free(u->image_data);
       exit(EXIT_FAILURE);
     }
   }


 }

void deallocate_image(image *u){

  int i;


  for(i = 0; i < u->m; i++){
    free(u->image_data[i]);
  }

  free(u->image_data);

}

void convert_jpeg_to_image(const unsigned char* image_chars, image *u){

    int i, j;


    for(i = 0; i < u->m ; i++){
      for(j = 0; j < u->n; j++){
        u->image_data[i][j] = (float)image_chars[i*u->n+j];
      }
    }

}


void convert_image_to_jpeg(const image *u, unsigned char* image_chars){

  int i,j;

  for(i = 0; i < u->m ; i++){
    for(j = 0; j < u->n; j++){
       image_chars[i*u->n+j]= (unsigned char)u->image_data[i][j];
    }
  }

}

/*
 * Velger å kjøre en iterasjon per iso_diffusion_denoising, slik at det blir lettere
 * med send og recv innad i prosessene
 */

void iso_diffusion_denoising(image *u, image *u_bar, float kappa){

  int i, j;
  int m = u->m;
  int n = u->n;


    for (i=1; i<m-1; i++)
        for (j=1; j<n-1; j++)
            u_bar->image_data[i][j] = u->image_data[i][j]
                + kappa*(u->image_data[i+1][j] + u->image_data[i][j+1]
                + u->image_data[i-1][j] + u->image_data[i][j-1]
                                         - 4*u->image_data[i][j]);

    for (i=1; i<m-1; i++)
        for (j=1; j<n-1; j++)
          u->image_data[i][j] = u_bar->image_data[i][j];

}
