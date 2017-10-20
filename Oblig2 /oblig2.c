#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <string.h>

typedef struct{

	double** matrix_elements;  // elements of matrix
	double* element_storage;   // contiguous storage
	int number_rows;           // y-axis
	int number_columns;	   // x-axis
} Mat;


void read_matrix_binaryformat(char* filename, Mat* matrix, int* num_rows, int* num_cols);
void write_matrix_binaryformat(char* filename, Mat matrix);
void allocate_matrix(Mat* matrix, int number_rows, int number_columns);
void deallocate_matrix(Mat* matrix);
void cannon_Multiply(int my_m, int my_l, int my_n, int my_rank,Mat matrix_a, Mat matrix_b, Mat matrix_c, MPI_Comm comm);
void multiply(Mat matrix_a, Mat matrix_b, Mat* matrix_c);
void print(Mat matrix);



int main(int argc, char* argv[]){

	int my_rank;
	int number_processes;
	int y_axisA; 						// Number of rows in A
	int x_axisA; 						// Number of columns in A
	int y_axisB; 						// Number of rows in B
	int x_axisB;  					// Number if columns in B

	int i, j, k; //counters used later

	Mat matrix_a;
	Mat matrix_b;
	Mat matrix_c;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
	MPI_Status status;

	int square;
	square = sqrt(number_processes - 1);   // Using this to check for right ammount of processes we need  in Canons algorithm
	int im_slave;
	if(my_rank != 0){
		im_slave = 1;
	}
	else{
		im_slave = MPI_UNDEFINED;
	}
	MPI_Comm slave_comm;
	MPI_Comm_split(MPI_COMM_WORLD, im_slave, my_rank, &slave_comm);

/*
 * Inside master process
 *
 * */
	if(my_rank == 0){
		char* infile_A;
		char* infile_B;
		char* outfile;

		Mat subA, subB, subC; // sub matrixes

		// Offsett values initialy 0
		int xoffsetA = 0;
		int yoffsetA = 0;
		int xoffsetB = 0;
		int yoffsetB = 0;


		if(argc != 4){
			fprintf(stderr, "Usage: arguments %s\n",argv[0]);
			exit(EXIT_FAILURE);
		}

		if(square*square != number_processes-1 ){
			fprintf(stderr, "Number of processes MUST must be square +1 %s\n",argv[0]);
			exit(EXIT_FAILURE);
		}

		infile_A = argv[1];
		infile_B = argv[2];
		outfile = argv[3];

		printf("outfile: %s\n", outfile);

		/* Lets read the binary files and allocate the right size of the outfile*/
		read_matrix_binaryformat(infile_A, &matrix_a, &y_axisA, &x_axisA);
		read_matrix_binaryformat(infile_B, &matrix_b, &y_axisB, &x_axisB);
		allocate_matrix(&matrix_c,y_axisA,x_axisB);
		printf("Allocated : %d x %d\n",y_axisA, x_axisB);


/*
* Testkode for å se om ting funker, det funker. Tallene jeg faar ut stemmer med egne beregninger
**/
		// y_axisA = 3;
		// x_axisA = 2;
		//
		// y_axisB = 2;
		// x_axisB = 3;
		//

		// allocate_matrix(&matrix_a, 3, 2);
		// allocate_matrix(&matrix_b, 2, 3);
		//
		// matrix_a.matrix_elements[0][0] = 1;
		// matrix_a.matrix_elements[0][1] = 2;
		// matrix_a.matrix_elements[1][0] = 3;
		// matrix_a.matrix_elements[1][1] = 4;
		// matrix_a.matrix_elements[2][0] = 5;
		// matrix_a.matrix_elements[2][1] = 6;
		//
		// matrix_b.matrix_elements[0][0] = 1;
		// matrix_b.matrix_elements[0][1] = 2;
		// matrix_b.matrix_elements[0][2] = 3;
		// matrix_b.matrix_elements[1][0] = 4;
		// matrix_b.matrix_elements[1][1] = 5;
		// matrix_b.matrix_elements[1][2] = 6;
		//
		// printf("Dette er A før noe:\n");
		// print(matrix_a);
		//
		// printf("Dette er B før noe:\n");
		//
		// print(matrix_b);

		//allocate_matrix(&matrix_c, 3, 3);

		//printf("%d vs %d\n", x_axisA, y_axisB);

		/* We  need to divide the matrixes into smaller parts and send it to the workers
		*  We also need to find the max sizes og the matrixes and if they are greater than the ammount of data we give them
		*  we set the reminnig values to zero
		*/
		int slave;

		for(slave = 1; slave < number_processes; slave ++){

			int Pi = ((slave -1) / square); // coordinates
			int Pj = ((slave -1) % square); // coordinates

			int my_m_A;      //  rows in A and
			int my_n_A;      //  columns in A that contains non zero values

			int my_m_B;     // rows in B and
			int my_n_B;	// columns in B tha contains non zero values
			int l_max;      // the inner two number in matrix (m,l) x (l,m)



			/*Rows in A*/
			if(Pi < (y_axisA%square)){
				my_m_A = y_axisA/square+1;
			}
			else{
				my_m_A = y_axisA/square;
			}
			/*Columns in A */

			if(Pj < (x_axisA%square)){
				my_n_A = x_axisA/square+1;
			}
			else{
				my_n_A = x_axisA/square;
			}

			/*Max size of 'l' in the whole matrix including zero values */

			if(0 < (x_axisA%square)){
				l_max = x_axisA/square+1;
			}
			else{
				l_max = x_axisA/square;
			}

			/*Rows in B*/

			if(Pi <(x_axisA%square)){
				my_m_B = x_axisA/square+1;
			}
			else{
				my_m_B = x_axisA/square;

			}

			/* Columns in B*/
			if(Pj <(x_axisB%square)){
				my_n_B = x_axisB/square+1;
			}
			else{
				my_n_B = x_axisB/square;

			}

			//printf("Fullsize %d vs realsize %d\n",l_max, my_n_A);

			/*Allocate the dub matrixes*/
			allocate_matrix(&subA, my_m_A, l_max);
			allocate_matrix(&subB, l_max, my_n_B);
			//printf("Just allocated matrixes\n subA :%d x %d\n and subB :%d x %d\n ",my_m_A, l_max, l_max, my_n_B);

			/*Filling in in the actuall size, the rest is left as zero values*/

			for(i = 0; i < my_m_A; i++){
				for(j = 0; j < my_n_A; j++){
					subA.matrix_elements[i][j] = matrix_a.matrix_elements[i+xoffsetA][j+yoffsetA];
				}
			}

			for(i = 0; i < my_m_B; i++){
				for(j = 0; j < my_n_B; j++){
					subB.matrix_elements[i][j] = matrix_b.matrix_elements[i+xoffsetB][j+yoffsetB];
				}
			}

			/* Oppdate the offsets. y-axis offsett are still 0 for any slave on far
 *right*/
			if(slave%square == 0){
				xoffsetA += my_m_A;
				xoffsetB += my_m_B;
				yoffsetA = 0;
				yoffsetB = 0;
			}
			else{
				yoffsetA +=my_n_A;
				yoffsetB +=my_n_B;
			}

			/*Now we are ready to send the sizes and subatrixes to the slaveren*/

			MPI_Send(&my_m_A, 1, MPI_INT, slave, 66, MPI_COMM_WORLD);
			MPI_Send(&l_max, 1, MPI_INT, slave, 66, MPI_COMM_WORLD);
			MPI_Send(&my_n_B, 1, MPI_INT, slave, 66, MPI_COMM_WORLD);

			MPI_Send(&subA.matrix_elements[0][0], (my_m_A*l_max),MPI_DOUBLE, slave,420,MPI_COMM_WORLD);
			MPI_Send(&subB.matrix_elements[0][0], (my_n_B*l_max) ,MPI_DOUBLE, slave,420,MPI_COMM_WORLD);

			/*Free the sub_matrixes*/
			deallocate_matrix(&subA);
			deallocate_matrix(&subB);

		} // end distrubute matrixes

		/*Clean up the large matrixes, we dont need them anymore*/

		deallocate_matrix(&matrix_a);
		deallocate_matrix(&matrix_b);



		/*Recieve the data from slave*/
		int xoffset_c = 0;
		int yoffset_c = 0;
		int my_m_C;
		int my_n_C;	// columns in C

		for(slave = 1; slave < number_processes; slave++){

			MPI_Recv(&my_m_C, 1, MPI_INT, slave, 666, MPI_COMM_WORLD, &status);
			MPI_Recv(&my_n_C, 1, MPI_INT, slave, 666, MPI_COMM_WORLD, &status);

			allocate_matrix(&subC, my_m_C, my_n_C);

			MPI_Recv(&subC.matrix_elements[0][0],my_m_C*my_n_C, MPI_DOUBLE, slave, 666, MPI_COMM_WORLD,&status);

			for(i = 0; i < my_m_C; i++){
				for(j = 0; j < my_n_C; j++){
					matrix_c.matrix_elements[i+xoffset_c][j+yoffset_c]= subC.matrix_elements[i][j];
				}

			}
			if(slave%square == 0 ){
				xoffset_c += my_m_C;
				yoffset_c = 0;
			}
			else{
				yoffset_c += my_n_C;
			}

			deallocate_matrix(&subC);

		}
		/*WRITE THE RESULT TO BINARY FILE*/

		//printf("Hello nr1\n");
		//print(matrix_c);
		write_matrix_binaryformat(outfile, matrix_c);
		deallocate_matrix(&matrix_c);

	} // end master

	/* Slave process
	 *
	 */
	else{
		int my_m, my_l, my_n;

		/*Recieveing the sizes of matix*/
		MPI_Recv(&my_m, 1, MPI_INT, 0, 66, MPI_COMM_WORLD,&status);
		MPI_Recv(&my_l, 1, MPI_INT, 0, 66, MPI_COMM_WORLD,&status);
		MPI_Recv(&my_n, 1, MPI_INT, 0, 66, MPI_COMM_WORLD,&status);

		/*Allocating matixes*/
		allocate_matrix(&matrix_a, my_m, my_l);
		allocate_matrix(&matrix_b, my_l, my_n);
		allocate_matrix(&matrix_c, my_m, my_n);


		/*Recieving the elements in matrixes*/
		MPI_Recv(&matrix_a.matrix_elements[0][0], (my_m*my_l), MPI_DOUBLE, 0, 420, MPI_COMM_WORLD,&status);
		MPI_Recv(&matrix_b.matrix_elements[0][0], (my_n*my_l), MPI_DOUBLE, 0, 420, MPI_COMM_WORLD,&status);
	//	printf("Hei fra main jeg er %d \n",my_rank);

	//	printf("Myrank Is %d\n",my_rank);



		/*Start multiplication. Credit goes to the text book */
		cannon_Multiply(my_m, my_l, my_n, my_rank, matrix_a, matrix_b,matrix_c, slave_comm);
	//	printf("Myrank Is %d\n",my_rank);


		/*Return result to master thread*/
		MPI_Send(&my_m, 1, MPI_INT, 0, 666, MPI_COMM_WORLD);
		MPI_Send(&my_n, 1, MPI_INT, 0, 666, MPI_COMM_WORLD);
		MPI_Send(&matrix_c.matrix_elements[0][0], (my_m*my_n), MPI_DOUBLE, 0, 666, MPI_COMM_WORLD);

//		MPI_Send(&matrix_c.matrix_elements[0][0],(my_m*my_n), MPI_DOUBLE, 0, 666, MPI_COMM_WORLD);

		/*The worker is done, lets dealloc*/
		deallocate_matrix(&matrix_a);
		deallocate_matrix(&matrix_b);
		deallocate_matrix(&matrix_c);


	}// Out of slave

	MPI_Finalize ();
	return 0;


}// End Main


void print(Mat matrix){

	int i;
	int j;

	for(i = 0; i < matrix.number_rows; i++){
		for(j=0; j < matrix.number_columns; j++){
			printf("[%d][%d] : %f\n",i,j,matrix.matrix_elements[i][j]);
		}
	}
	for(i = 0; i < matrix.number_rows*matrix.number_columns; i++){

			printf("[%d] : %f\n",i,matrix.element_storage[i]);

	}
}



/*
 * Taken straigth from the textbook full credit to the book there.
 * **/

void cannon_Multiply(int my_m, int my_l, int my_n, int my_rank, Mat matrix_a, Mat matrix_b, Mat matrix_c, MPI_Comm comm){

	int i, j;
	int npes;
	int dims[2];
	int periods[2];
	int myrank;
	int my2drank;
	int mycoords[2];
	int uprank;
	int downrank;
	int leftrank;
	int rightrank;
	int coords[2];
	int shiftsource;
	int shiftdest;
	MPI_Status status;
	MPI_Comm comm_2d;

	/*Get the communicator related info*/
	MPI_Comm_size(comm, &npes);
	MPI_Comm_rank(comm, &myrank);

//	printf("Vi er egentlig %d stykker\n",npes);

//	printf("Forste test om %d er her \n",myrank);

	/*Set up the Cartesian topology*/

	dims[0] = dims[1] = sqrt(npes);

	/*Set up the periods for wraparound connections */
	periods[0] = periods[1] = 1;


	/*Create the Catesian topology, with rank reordering*/
	MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);


	/*Get the rank and coordinates with repect to new topology*/
	MPI_Comm_rank(comm_2d, &my2drank);
	MPI_Cart_coords(comm_2d, my2drank, 2 ,mycoords);

	/*Compute ranks of the up and left shifts*/
	MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
	MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

	/*Perform the initial matrix alignment. First for A then for B*/
	MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(&matrix_a.matrix_elements[0][0], my_m*my_l, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(&matrix_b.matrix_elements[0][0], my_n*my_l, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	/*Start the main computation loop*/

//	printf("%d kom hit\n",myrank);

	for(i = 0; i < dims[0]; i++){
		multiply(matrix_a, matrix_b, &matrix_c);

		/*shift matrix left by one*/
		MPI_Sendrecv_replace(&matrix_a.matrix_elements[0][0], my_m*my_l, MPI_DOUBLE, leftrank, 1, rightrank, 1, comm_2d, &status);

		/*shift matrix b up by one*/

		MPI_Sendrecv_replace(&matrix_b.matrix_elements[0][0], my_n*my_l, MPI_DOUBLE, uprank, 1, downrank, 1, comm_2d, &status);

	}

	/* Free the communicator*/

	MPI_Comm_free(&comm_2d);


}

/*Multiplication method that are done using OpenMP*/

void multiply(Mat matrix_a, Mat matrix_b, Mat* matrix_c){
	int i,j,k;

	int m = matrix_a.number_rows;
	int l = matrix_a.number_columns;
	int n = matrix_b.number_columns;

	#pragma omp parallel for private(j,k)
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k <l;k++){
			matrix_c->matrix_elements[i][j] += matrix_a.matrix_elements[i][k]*matrix_b.matrix_elements[k][j];
			}
		}
	}
//	printf("HELOOOO\n");
}


/* Reads a matrix in binary form
 * and stores it in a matrix struct
 * makes call to func: allocate_matrix(Mat* matrix, int number_rows, int number_columns)
 * */


void read_matrix_binaryformat(char* filename, Mat* matrix, int* num_rows, int* num_columns){

	FILE* fp = fopen(filename,"rb");
	fread(num_rows, sizeof(int), 1, fp);
	fread(num_columns, sizeof(int), 1, fp);

	/*storage allocation for the matrix*/
	allocate_matrix(matrix, *num_rows, *num_columns);

	/*read in the entire matrix*/
	fread((*matrix).matrix_elements[0],sizeof(double), (*num_rows)*(*num_columns),fp);
	close(fp);
}


/*
 * Writes result of matrix multiplication to a file in binary form.
 * m and n are really just the same integer since we do matrix-multiplication
 * */

void write_matrix_binaryformat(char* filename, Mat matrix){

	int m = matrix.number_rows;
	int n = matrix.number_columns;
	printf("Hello\n");

	FILE* fp = fopen(filename,"wb");
	fwrite(&(matrix.number_rows), sizeof(int), 1, fp);
	fwrite(&(matrix.number_columns), sizeof(int), 1, fp);
	fwrite(matrix.matrix_elements[0], sizeof(double),m*n ,fp);

	fclose(fp);
}

/**
 * Allocates memory for a matrix struct.
 * Note to self: should I set all elements to 0 in an inner for-loop?
 * */


void allocate_matrix(Mat* matrix, int num_rows, int num_columns){

	int i;
	matrix->number_rows = num_rows;
	matrix->number_columns =  num_columns;

	if((matrix->matrix_elements = (double**)malloc(num_rows*sizeof(double*)))== NULL){
		perror("Failed to malloc matrix_elements");
		exit(EXIT_FAILURE);
	}

	if((matrix->element_storage = (double*)malloc(num_rows*num_columns*sizeof(double)))==NULL){
		perror("Failed to malloc element_storage");
		free(matrix->matrix_elements);
		exit(EXIT_FAILURE);
	}

	for(i = 0; i <num_rows; i++){
		matrix->matrix_elements[i] = &matrix->element_storage[i*num_columns];
	}

}

void deallocate_matrix(Mat* matrix){

	free(matrix->matrix_elements);
	free(matrix->element_storage);
}
