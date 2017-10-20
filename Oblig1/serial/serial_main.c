#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {

  float** image_data;
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
void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters);



int main(int argc, char *argv[]) {

    /*
      Decalration of variables
    */

    int m, n, c, iters;
    float kappa;
    image u, u_bar;
    unsigned char *image_chars;
    char *input_jpeg_filename;
    char *output_jpeg_filename;


    /*
      Checking if input is correct
    */

    if(argc != 5){
      fprintf(stderr, "Usage: arguments %s\n",argv[0]);
      exit(EXIT_FAILURE);
    }

    if((kappa = atof(argv[1])) <= 0){
      fprintf(stderr, "Usage: kappa too small or invalid : '%s'\n",argv[1]);
      exit(EXIT_FAILURE);
    }
    printf(" kappa : %f \n", kappa);


    if((iters = atoi(argv[2])) <= 0){
      fprintf(stderr, "Usage: iter too small or invalid : '%s'\n",argv[2]);
      exit(EXIT_FAILURE);
    }

    printf(" Iters : %d\n",iters);

    input_jpeg_filename = argv[3];
    printf("input name : %s\n", input_jpeg_filename);

    output_jpeg_filename = argv[4];
    printf("output name : %s\n", output_jpeg_filename);


    import_JPEG_file(input_jpeg_filename, &image_chars,&m ,&n ,&c);

    allocate_image (&u, m, n);

    allocate_image(&u_bar, m, n);

    convert_jpeg_to_image (image_chars, &u);


    iso_diffusion_denoising (&u, &u_bar, kappa, iters);


    convert_image_to_jpeg (&u_bar, image_chars);


    export_JPEG_file(output_jpeg_filename, image_chars, m, n, 1, 75);
    deallocate_image(&u);
    deallocate_image(&u_bar);





  return 0;
}
/*
 * @param
 *
 *
 * Exits on any failure on systemscalls
*/

 void allocate_image(image *u, int m, int n){


   memcpy(&u->m,&m,sizeof(int));
   memcpy(&u->n,&n,sizeof(int));

//   printf("Memcpy m : %d\n",m );
  //  u-> m = m;
  //  u-> n = n;

   printf("Size on image: %d x %d\n",u->m, u->n );


   if((u->image_data = (float**)malloc(m * sizeof(float*))) == NULL){
    perror("Failed malloc pointer to matrix");
    exit(EXIT_FAILURE);
   }

   int i;
   //Test om ending skjer

   for(i = 0; i < m; i++){
     if((u->image_data[i] = (float*)malloc(n * sizeof(float))) == NULL){
       perror("Inner malloc");
       free(u->image_data);
       exit(EXIT_FAILURE);
     }
   }

   printf("Test of might : %d x %d \n",u->m , u->n );

 }

void deallocate_image(image *u){

  int i;

  printf("My m in deallocate_image : %d\n",u->m);
  printf("My n in deallocate_image : %d \n",u->n );

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

void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters){

  int i, j, counter;
  counter = 0;
  printf("Hello from iso \n" );
  printf("u_bar->image_data[0][0] : %f\n",u_bar->image_data[4000][2000]);
  printf("u->image_data[0][0] : %f\n",u->image_data[4000][2000]);
  printf("u->m :%d\n",u->m);
  printf("u->n :%d\n",u->n);


  while(counter < iters ){

    for(i = 0; i < u->m ; i++){
      for(j = 0; j < u->n ; j++){
        // printf("Inner for\n" );
        // printf("i : %d J : %d \n",i,j );
        if((1 <= i) &&  (i <= (u->m -2)) && (1 <= j) && (j <= (u->n-2))){
          // printf("if\n" );
          u_bar->image_data[i][j] = u->image_data[i][j] + kappa*(u->image_data[i-1][j] +
             u->image_data[i][j-1] + u->image_data[i][j+1] + u->image_data[i+1][j] -
              4*u->image_data[i][j]);
        }
        else{
          // printf("else\n");
          u_bar->image_data[i][j] = u->image_data[i][j]; //Bounday pixels are copied into u_bar
        }
      }
    }

    //TODO: se om det er mulig å kopiere tilbake i samme for-løkke som over

      for(i = 0 ; i < u->m ; i++){
        for(j = 0; j < u->n ; j++){
          u->image_data[i][j] = u_bar->image_data[i][j];
      }
    }

    counter ++;
  } // end while

}
