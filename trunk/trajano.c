/*/
 * Axel S. Trajano
 * 2002-00283
 * CMSC 190 (2) : Special Problem
 * Implementation and Evaluation of Selected Image Processing Routines
 * for Single and Distributed Processors
 * $Id$
/*/

#include <math.h>
#include <stdio.h>
#include <jpeglib.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <mpi.h>



int numprocs, rank;
int namelen;
/* we will be using this uninitialized pointer later to store raw, uncompressd image */
unsigned char *raw_image = NULL;


/* dimensions of the image we want to write */
int width;
int height;
int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */
int cspace, icomp;
//char st[];

int rem, quo, start, end;

#define rout 12
/* Time monitoring purposes */
double Totalroutine[rout];	// total time for a specific routine
double Avgroutine[rout];	// average time of the execution of a specific routine
double Exectime = 0;		// Overall execution time
double AvgExec = 0;			// Average execution time
int rcnt[rout];

char routines[12][30] = {"Mean","Median","Brightness","Contrast","Invert",
			"Threshold","Noise Reduction","Laplace (4-N)",
			"Laplace (8-N)","Sobel","Dilation","Erosion"};


/**
 * read_jpeg_file Reads from a jpeg file on disk specified by filename and saves into the 
 * raw_image buffer in an uncompressed format.
 * 
 * \returns positive integer if successful, -1 otherwise
 * \param *filename char string specifying the file name to read from
 *
 */

int read_jpeg_file( char *filename )
{


	/* these are standard libjpeg structures for reading(decompression) */
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	/* libjpeg data structure for storing one row, that is, scanline of an image */
	JSAMPROW row_pointer[1];
	
	FILE *infile = fopen( filename, "rb" );
	unsigned long location = 0;
	int i = 0;
	

	for(i=0;i<rout;i++){
		Totalroutine[i] = 0;
		Avgroutine[i] = 0;
		rcnt[i] = 0;
	}

	if ( !infile )
	{
		printf("Error opening jpeg file %s! \n", filename );
		return -1;
	}

	/* here we set up the standard libjpeg error handler */
	cinfo.err = jpeg_std_error( &jerr );
	/* setup decompression process and source, then read JPEG header */
	jpeg_create_decompress( &cinfo );
	/* this makes the library read from infile */
	jpeg_stdio_src( &cinfo, infile );
	/* reading the image header which contains image information */
	jpeg_read_header( &cinfo, TRUE );

	/* Output image information, if needed. */

   if(rank == 0){
	printf("\n");	
	printf( "JPEG File Information: \n" );
	printf( "Image width and height: %d pixels and %d pixels.\n", cinfo.image_width, cinfo.image_height );
	printf( "Color components per pixel: %d.\n", cinfo.num_components );
	printf( "Color space: %d.\n\n", cinfo.jpeg_color_space );
	
   }
	icomp = cinfo.num_components;
	width = cinfo.image_width;
	height = cinfo.image_height;
	unsigned char *samp;	

	int j;
	


	/* This will force the read color space to be grayscaled */
	cinfo.out_color_space = JCS_GRAYSCALE;
	
	
	cspace = cinfo.out_color_space;
	/* Start decompression jpeg here */
	jpeg_start_decompress( &cinfo );

	/* allocate memory to hold the uncompressed image */
	raw_image = (unsigned char*)malloc( cinfo.output_width*cinfo.output_height*cinfo.num_components );
	/* now actually read the jpeg into the raw buffer */
	row_pointer[0] = (unsigned char *)malloc( cinfo.output_width*cinfo.num_components );


	samp = (unsigned char *)malloc(cinfo.output_height*cinfo.output_width*cinfo.num_components*sizeof(unsigned char));

	int loc = 0;
	

	/* read one scan line at a time */
	while( cinfo.output_scanline < cinfo.image_height )
	{
		jpeg_read_scanlines( &cinfo, row_pointer, 1 );
		for( i=0; i<cinfo.image_width*cinfo.num_components;i++){ 
			//matrix[loc][i] = row_pointer[0][i];
			samp[loc*cinfo.image_width+i] = row_pointer[0][i];
			raw_image[location++] = row_pointer[0][i];
		}
		loc++;
		
	}

	loc = 0;

		
	int swtch = 100;
	int ok;
	char opt[800];
	
	strcpy(opt,"");


	do{

	   if(rank == 0){
		printf("Select Operators: \n");
		printf("\t 1. Mean \n");
		printf("\t 2. Median \n");
		printf("\t 3. Brightness \n");
		printf("\t 4. Contrast \n");
		printf("\t 5. Invert \n");
		printf("\t 6. Threshold \n");
		printf("\t 7. Noise Reduction \n");
		printf("\t 8. Laplace (4-Neighbor) \n");
		printf("\t 9. Laplace (8-Neighbor) \n");
		printf("\t 10. Sobel \n");
		printf("\t 11. Dilation \n");
		printf("\t 12. Erosion \n");
		printf("\t 0. Exit \n");
		
		printf("Choice:\t");
		scanf("%d", &swtch);


	   } // end of root access


		MPI_Bcast(&swtch,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		

		switch(swtch){

			case 1: mean_seq(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[mean] "); rcnt[0]++;
				Avgroutine[0] = Totalroutine[0] / rcnt[0];
				break;
			case 2: get_median(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[median] "); rcnt[1]++;
				Avgroutine[1] = Totalroutine[1] / rcnt[1];
				break;
			case 3: brightness(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[brightness] "); rcnt[2]++;
				Avgroutine[2] = Totalroutine[2] / rcnt[2];
				break;
			case 4: contrast(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[contrast] "); rcnt[3]++;
				Avgroutine[3] = Totalroutine[3] / rcnt[3];
				break;
			case 5: invert(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[invert] "); rcnt[4]++;
				Avgroutine[4] = Totalroutine[4] / rcnt[4];
				break;
			case 6: threshold(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[threshold] "); rcnt[5]++;
				Avgroutine[5] = Totalroutine[5] / rcnt[5];
				break;
			case 7: noise_reduce(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[noise reduce] "); rcnt[6]++;
				Avgroutine[6] = Totalroutine[6] / rcnt[6];
				break;
			case 8: edge_laplace(samp, cinfo.image_width, cinfo.num_components, 0);
				strcat(opt,"[(4-neighbor) laplace] "); rcnt[7]++;
				Avgroutine[7] = Totalroutine[7] / rcnt[7];
				break;
			case 9: edge_laplace(samp, cinfo.image_width, cinfo.num_components, 1);
				strcat(opt,"[(8-neighbor) laplace] "); rcnt[8]++;
				Avgroutine[8] = Totalroutine[8] / rcnt[8];
				break;
			case 10: edge_sobel(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[sobel] "); rcnt[9]++;
				Avgroutine[9] = Totalroutine[9] / rcnt[9];
				break;
			case 11: dilation(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[dilate] "); rcnt[10]++;
				Avgroutine[10] = Totalroutine[10] / rcnt[10];
				break;
			case 12: erosion(samp, cinfo.image_width, cinfo.num_components);
				strcat(opt,"[erode] "); rcnt[11]++;
				Avgroutine[11] = Totalroutine[11] / rcnt[11];
				break;
			case 0: 
			default: break;

		}
		//system("clear");

	}while(swtch != 0);

   if(rank == 0)
	printf("\nOperations: \n%s \n\n",opt);
	

	/* wrap up decompression, destroy objects, free pointers and close open files */
	jpeg_finish_decompress( &cinfo );
	jpeg_destroy_decompress( &cinfo );
	free( row_pointer[0] );
	

	free(samp);


	fclose( infile );
	/* yup, we succeeded! */
	return 1;

}


/*=======================================================================================*/
/*======================================== MEDIAN =======================================*/
/*=======================================================================================*/

unsigned char median(unsigned char a, unsigned char b, unsigned char c){

	unsigned char temp;

	if(a > b){
		temp = a; a = b; b = temp;
	}
	if(b > c){
		temp = b; b = c; c = temp;
	}
	if(a > b){
		temp = a; a = b; b = temp;
	}

	return b;

}


int get_median(unsigned char *matrix, int w, int comp){

	double time;
	int row = 0;
	int i, j, x, y, shift;
	unsigned long r, c, rsize, rsize2, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], disp2[numprocs], nsent[numprocs], nsent2[numprocs];
	unsigned long location = 0;
	unsigned char mid, aa, bb, cc;

	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;


	if((rank == 0) || (rank == numprocs-1))
		rsize2 = rsize + 1; // rank 0 has lower row, and last rank has upper row;
	else
		rsize2 = rsize + 2; // middle ranks has upper and lower rows;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));
	unsigned char *submat = (unsigned char *)malloc(rsize2*w*comp*sizeof(unsigned char));


	// set the buffer size
	unsigned long bsize = rsize * w;
	unsigned long bsize2 = rsize2 * w;


	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}
		// Values for submat[][]
		for(i=0;i<numprocs;i++){
			if(i!=0)
				disp2[i] = disp[i] - w;
			else
				disp2[i] = disp[i];
			if(i == numprocs-1)
				nsent2[i] = nsent[i] + w;
			else{
				nsent2[i] = nsent[i] + (2 * w);
				if(i==0)
					nsent2[i] = nsent[i] + w;
			}
		}
		

	} // end of root task


	if(rank != 0)
		row = 1;

	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	MPI_Scatterv(matrix,nsent2,disp2,MPI_UNSIGNED_CHAR,submat,bsize2,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);	
	
	//printf("scattered \n");

	time = MPI_Wtime();

	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			
			if( ((rank==0) && (r==0)) || ((rank==numprocs-1) && (r==rsize)) || c==0 || c==width-1 )
				mid = 122;
				
			else{
				aa = median(submat[((r+row-1)*w)+(c-1)],submat[((r+row-1)*w)+c],submat[((r+row-1)*w)+(c+1)]);
				bb = median(submat[((r+row)*w)+(c-1)],submat[((r+row)*w)+c],submat[((r+row)*w)+(c+1)]);
				cc = median(submat[((r+row+1)*w)+(c-1)],submat[((r+row+1)*w)+c],submat[((r+row+1)*w)+(c+1)]);

				mid = median(aa, bb, cc);

				A[(r*w)+c] = mid;
			}

		} // end of column traversal
	} // end of row traversal

	time = MPI_Wtime() - time;
	
	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);


	Totalroutine[1] += time;

	//printf("Gathered rank %d \n",rank);


	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A); free(submat);

	return 1;

}


/*=======================================================================================*/
/*===================================== BRIGHTNESS ======================================*/
/*=======================================================================================*/

int brightness(unsigned char *matrix, int w, int comp){

	double time;
	int row = 0;
	int i, j, x, y, shift;
	unsigned long r, c, rsize, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], nsent[numprocs];
	unsigned long location = 0;
	unsigned char br;
	float bright;

	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));

	// set the buffer size
	unsigned long bsize = rsize * w;;

	if(rank == 0){
		printf("Enter a numeric Brightness value: ");
		scanf("%f",&bright);
	}

	MPI_Bcast(&bright,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}

	} // end of root task


	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	//printf("scattered \n");

	time = MPI_Wtime();

	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			br = A[(r*w)+c] + bright;
			if(br < 0)
				br = 0;
			else if(br > 255)
				br = 255;
			A[(r*w)+c] = br;
		}
	}

	time = MPI_Wtime() - time;

	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);


	Totalroutine[2] += time;



	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A);


	return 1;

}

/*=======================================================================================*/
/*======================================= CONTRAST ======================================*/
/*=======================================================================================*/

int contrast(unsigned char *matrix, int w, int comp){

	double time;
	int row = 0;
	int i, j, x, y, shift;
	unsigned long r, c, rsize, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], nsent[numprocs];
	unsigned long location = 0;
	unsigned char ct;
	float contrast;

	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));

	// set the buffer size
	unsigned long bsize = rsize * w;

	if(rank == 0){
		printf("Enter a numeric Contrast value: ");
		scanf("%f",&contrast);
	}

	MPI_Bcast(&contrast,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}

	} // end of root task


	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	//printf("scattered \n");

	time = MPI_Wtime();

	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			ct = ((A[(r*w)+c]-0.5) * contrast) + 0.5;
			if(ct < 0)
				ct = 0;
			else if(ct > 255)
				ct = 255;
			A[(r*w)+c] = ct;
		}
	}

	time = MPI_Wtime() - time;

	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);


	Totalroutine[3] += time;



	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A);


	return 1;

}


/*=======================================================================================*/
/*========================================= INVERT ======================================*/
/*=======================================================================================*/

int invert(unsigned char *matrix, int w, int comp){

	double time;
	int row = 0;
	int i, j, x, y, shift;
	unsigned long r, c, rsize, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], nsent[numprocs];
	unsigned long location = 0;
	int inv;

	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));


	// set the buffer size
	unsigned long bsize = rsize * w;


	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}	

	} // end of root task


	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	//printf("scattered \n");

	time = MPI_Wtime();

	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){
			inv = 1 - A[(r*w)+c];
			A[(r*w)+c] = inv;
		}
	}

	time = MPI_Wtime() - time;

	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);


	Totalroutine[4] += time;


	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A);



	return 1;

}



/*=======================================================================================*/
/*===================================== THRESHOLD =======================================*/
/*=======================================================================================*/

int threshold(unsigned char *matrix, int w, int comp){

	double time;
	int row = 0;
	int i, j, x, y, shift;
	unsigned long r, c, rsize, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], nsent[numprocs];
	unsigned long location = 0;
	unsigned int thresh;

	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));


	// set the buffer size
	unsigned long bsize = rsize * w;

	if(rank == 0){
		do{
			printf("Enter a Threshold value (0 - 255): ");
			scanf("%d",&thresh);
		}while(thresh < 0 || thresh > 255);
	}

	MPI_Bcast(&thresh,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);


	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}

	} // end of root task


	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	//printf("scattered \n");

	time = MPI_Wtime();

	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			if(A[(r*w)+c] < thresh)
				A[(r*w)+c] = 0;
			else
				A[(r*w)+c] = 255;
		}
	}

	time = MPI_Wtime() - time;

	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);


	Totalroutine[5] += time;



	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A);

	return 1;

}


/*=======================================================================================*/
/*================================= DILATION & EROSION ==================================*/
/*=======================================================================================*/

int dilation(unsigned char *matrix, int w, int comp){


	double time;
	int row = 0;
	int i, j, x, y, shift;
	unsigned long r, c, rsize, rsize2, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], disp2[numprocs], nsent[numprocs], nsent2[numprocs];
	unsigned long location = 0;
	boolean check = 0;

	threshold(matrix, w, comp);

	int mask[3][3];

	// 0 represents the black's
	mask[0][0] = 0; mask[0][1] = 0; mask[0][2] = 0;
	mask[1][0] = 0; mask[1][1] = 0; mask[1][2] = 0;
	mask[2][0] = 0; mask[2][1] = 0; mask[2][2] = 0;


	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;


	if((rank == 0) || (rank == numprocs-1))
		rsize2 = rsize + 1; // rank 0 has lower row, and last rank has upper row;
	else
		rsize2 = rsize + 2; // middle ranks has upper and lower rows;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));
	unsigned char *submat = (unsigned char *)malloc(rsize2*w*comp*sizeof(unsigned char));


	// set the buffer size
	unsigned long bsize = rsize * w;
	unsigned long bsize2 = rsize2 * w;


	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}
		// Values for submat[][]
		for(i=0;i<numprocs;i++){
			if(i!=0)
				disp2[i] = disp[i] - w;
			else
				disp2[i] = disp[i];
			if(i == numprocs-1)
				nsent2[i] = nsent[i] + w;
			else{
				nsent2[i] = nsent[i] + (2 * w);
				if(i==0)
					nsent2[i] = nsent[i] + w;
			}
		}
		

	} // end of root task

	if(rank != 0)
		row = 1;

	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	MPI_Scatterv(matrix,nsent2,disp2,MPI_UNSIGNED_CHAR,submat,bsize2,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);	
	
	//printf("scattered \n");

	time = MPI_Wtime();

	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			int sum = 0;
			if( ((rank==0) && (r==0)) || ((rank==numprocs-1) && (r==rsize)) || c==0 || c==width-1 )
				sum = 0;
				
			else{

				check = 0;
				for(x=-1; x<2; x++){
					for(y=-1; y<2; y++){
						if(mask[x+1][y+1] == 0)
						if(submat[((r+row+x)*w)+(c+y)] == mask[x+1][y+1]){
							check = 1;
							break;
						}
					}
					if(check)
						break;
					
				}
				if(check){
					A[(r*w)+c] = 0;
				}
				else
					A[(r*w)+c] = 255;
				
				check = 0;

			}

		} // end of column traversal

	} // end of row traversal


	for(r=0; r<rsize; r++)
		for(c=0; c<width; c++)
			if(A[(r*w)+c]==0 || submat[((r+row)*w)+c]==0)
				A[(r*w)+c] = 0;


	time = MPI_Wtime() - time;

	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);


	Totalroutine[10] += time;




	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A); free(submat);


	return 1;

}



int erosion(unsigned char *matrix, int w, int comp){

	double time;
	int row = 0;
	int i, j, x, y, shift;
	unsigned long r, c, rsize, rsize2, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], disp2[numprocs], nsent[numprocs], nsent2[numprocs];
	unsigned long location = 0;
	boolean check = 0;

	threshold(matrix, w, comp);

	int mask[3][3];

	// 0 represents the black's
	mask[0][0] = 0; mask[0][1] = 0; mask[0][2] = 0;
	mask[1][0] = 0; mask[1][1] = 0; mask[1][2] = 0;
	mask[2][0] = 0; mask[2][1] = 0; mask[2][2] = 0;


	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;


	if((rank == 0) || (rank == numprocs-1))
		rsize2 = rsize + 1; // rank 0 has lower row, and last rank has upper row;
	else
		rsize2 = rsize + 2; // middle ranks has upper and lower rows;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));
	unsigned char *submat = (unsigned char *)malloc(rsize2*w*comp*sizeof(unsigned char));


	// set the buffer size
	unsigned long bsize = rsize * w;
	unsigned long bsize2 = rsize2 * w;


	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}
		// Values for submat[][]
		for(i=0;i<numprocs;i++){
			if(i!=0)
				disp2[i] = disp[i] - w;
			else
				disp2[i] = disp[i];
			if(i == numprocs-1)
				nsent2[i] = nsent[i] + w;
			else{
				nsent2[i] = nsent[i] + (2 * w);
				if(i==0)
					nsent2[i] = nsent[i] + w;
			}
		}
		

	} // end of root task

	if(rank != 0)
		row = 1;

	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	MPI_Scatterv(matrix,nsent2,disp2,MPI_UNSIGNED_CHAR,submat,bsize2,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);	
	
	//printf("scattered \n");

	time = MPI_Wtime();

	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			int sum = 0;
			if( ((rank==0) && (r==0)) || ((rank==numprocs-1) && (r==rsize)) || c==0 || c==width-1 )
				sum = 0;
				
			else{

				check = 0;
				for(x=-1; x<2; x++){
					for(y=-1; y<2; y++){
						if(mask[x+1][y+1] == 0)
						if(submat[((r+row+x)*w)+(c+y)] == mask[x+1][y+1]){
							check = 1;
							break;
						}
					}
					if(check)
						break;
					
				}
				if(check){
					A[(r*w)+c] = 0;
				}
				else
					A[(r*w)+c] = 255;
				
				check = 0;

			}

		} // end of column traversal

	} // end of row traversal

	for(r=0; r<rsize; r++)
		for(c=0; c<width; c++)
			if(A[(r*w)+c]==0 && submat[((r+row)*w)+c]==0)
				A[(r*w)+c] = 0;

	time = MPI_Wtime() - time;

	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);



	Totalroutine[11] += time;




	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A); free(submat);


	return 1;


}


/*=======================================================================================*/
/*================================== NOISE REDUCTION ====================================*/
/*=======================================================================================*/

int noise_reduce(unsigned char *matrix, int w, int comp){


	double time;
	int row = 0;
	int i, j, x, y, shift;
	unsigned long r, c, rsize, rsize2, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], disp2[numprocs], nsent[numprocs], nsent2[numprocs];
	unsigned long location = 0;


	int mask[3][3];

	
	mask[0][0] = 1; mask[0][1] = 1; mask[0][2] = 1;
	mask[1][0] = 1; mask[1][1] = 8; mask[1][2] = 1;
	mask[2][0] = 1; mask[2][1] = 1; mask[2][2] = 1;


	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;


	if((rank == 0) || (rank == numprocs-1))
		rsize2 = rsize + 1; // rank 0 has lower row, and last rank has upper row;
	else
		rsize2 = rsize + 2; // middle ranks has upper and lower rows;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));
	unsigned char *submat = (unsigned char *)malloc(rsize2*w*comp*sizeof(unsigned char));


	// set the buffer size
	unsigned long bsize = rsize * w;
	unsigned long bsize2 = rsize2 * w;



	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}
		// Values for submat[][]
		for(i=0;i<numprocs;i++){
			if(i!=0)
				disp2[i] = disp[i] - w;
			else
				disp2[i] = disp[i];
			if(i == numprocs-1)
				nsent2[i] = nsent[i] + w;
			else{
				nsent2[i] = nsent[i] + (2 * w);
				if(i==0)
					nsent2[i] = nsent[i] + w;
			}
		}
		

	} // end of root task

	if(rank != 0)
		row = 1;

	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	MPI_Scatterv(matrix,nsent2,disp2,MPI_UNSIGNED_CHAR,submat,bsize2,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);	
	
	//printf("scattered \n");


	time = MPI_Wtime();

	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			int sum = 0;
			if( ((rank==0) && (r==0)) || ((rank==numprocs-1) && (r==rsize)) || c==0 || c==width-1 )
				sum = 0;
				
			else{
				for(x=-1; x<2; x++){
					for(y=-1; y<2; y++){
						sum = sum + submat[((r+row+x)*w)+(c+y)] * mask[x+1][y+1];
					}
				}
				
				sum = sum / 9;
				if(sum < 0)
					sum = 0;
				else if(sum > 255)
					sum = 255;
				A[(r*w)+c] = sum;

			}

		} // end of column traversal
	} // end of row traversal

	time = MPI_Wtime() - time;

	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);



	Totalroutine[6] += time;




	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A); free(submat);


	return 1;

}


/*=======================================================================================*/
/*======================================== MEAN =========================================*/
/*=======================================================================================*/

int mean_seq(unsigned char *matrix, int w, int comp){

	double time;
	int row = 0;
	int i, j, x, y, shift;
	unsigned long r, c, rsize, rsize2, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], disp2[numprocs], nsent[numprocs], nsent2[numprocs];
	unsigned long location = 0;

/*
	int mask[3][3];

	// can be done without this mask.
	mask[0][0] = 1; mask[0][1] = 1; mask[0][2] = 1;
	mask[1][0] = 1; mask[1][1] = 1; mask[1][2] = 1;
	mask[2][0] = 1; mask[2][1] = 1; mask[2][2] = 1;
*/

	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;


	if((rank == 0) || (rank == numprocs-1))
		rsize2 = rsize + 1; // rank 0 has lower row, and last rank has upper row;
	else
		rsize2 = rsize + 2; // middle ranks has upper and lower rows;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));
	unsigned char *submat = (unsigned char *)malloc(rsize2*w*comp*sizeof(unsigned char));


	// set the buffer size
	unsigned long bsize = rsize * w;
	unsigned long bsize2 = rsize2 * w;



	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}
		// Values for submat[][]
		for(i=0;i<numprocs;i++){
			if(i!=0)
				disp2[i] = disp[i] - w;
			else
				disp2[i] = disp[i];
			if(i == numprocs-1)
				nsent2[i] = nsent[i] + w;
			else{
				nsent2[i] = nsent[i] + (2 * w);
				if(i==0)
					nsent2[i] = nsent[i] + w;
			}
		}
		

	} // end of root task

	if(rank != 0)
		row = 1;

	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	MPI_Scatterv(matrix,nsent2,disp2,MPI_UNSIGNED_CHAR,submat,bsize2,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);	
	
	//printf("scattered \n");

	time = MPI_Wtime();

	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			int sum = 0;
			if( ((rank==0) && (r==0)) || ((rank==numprocs-1) && (r==rsize)) || c==0 || c==width-1 )
				sum = 0;
				
			else{
				for(x=-1; x<2; x++){
					for(y=-1; y<2; y++){
						sum = sum + submat[((r+row+x)*w)+(c+y)]; //* mask[x+1][y+1];
					}
				}
				
				sum = sum / 9;
				if(sum < 0)
					sum = 0;
				else if(sum > 255)
					sum = 255;

				A[(r*w)+c] =  sum;

			}

		} // end of column traversal

	} // end of row traversal
	time = MPI_Wtime() - time;


	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	//printf("Gathered rank %d \n",rank);


	Totalroutine[0] += time;



	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A); free(submat);

	return 1;

}


/*=======================================================================================*/
/*======================================= LAPLACE =======================================*/
/*=======================================================================================*/

int edge_laplace(unsigned char *matrix, int w, int comp, int n){


	double time;
	int row = 0;
	int i, j, x, y, index, shift;
	unsigned long r, c, rsize, rsize2, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], disp2[numprocs], nsent[numprocs], nsent2[numprocs];
	unsigned long location = 0;


	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;

	int mask[3][3];

if(n == 0){
	// 4-Neighbor
	mask[0][0] = 0; mask[0][1] =-1; mask[0][2] = 0;
	mask[1][0] =-1; mask[1][1] = 4; mask[1][2] =-1;
	mask[2][0] = 0; mask[2][1] =-1; mask[2][2] = 0;
	index = 7;
}	
else{
	// 8-Neighbor
	mask[0][0] =-1; mask[0][1] =-1; mask[0][2] =-1;
	mask[1][0] =-1; mask[1][1] = 8; mask[1][2] =-1;
	mask[2][0] =-1; mask[2][1] =-1; mask[2][2] =-1;
	index = 8;
}

	if((rank == 0) || (rank == numprocs-1))
		rsize2 = rsize + 1; // rank 0 has lower row, and last rank has upper row;
	else
		rsize2 = rsize + 2; // middle ranks has upper and lower rows;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));
	unsigned char *submat = (unsigned char *)malloc(rsize2*w*comp*sizeof(unsigned char));


	// set the buffer size
	unsigned long bsize = rsize * w;
	unsigned long bsize2 = rsize2 * w;



	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}
		// Values for submat[][]
		for(i=0;i<numprocs;i++){
			if(i!=0)
				disp2[i] = disp[i] - w;
			else
				disp2[i] = disp[i];
			if(i == numprocs-1)
				nsent2[i] = nsent[i] + w;
			else{
				nsent2[i] = nsent[i] + (2 * w);
				if(i==0)
					nsent2[i] = nsent[i] + w;
			}
		}
		

	} // end of root task

	
	if(rank != 0)
		row = 1;


	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	MPI_Scatterv(matrix,nsent2,disp2,MPI_UNSIGNED_CHAR,submat,bsize2,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);	
	
	//printf("scattered \n");


	time = MPI_Wtime();


	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			int sum = 0;
			if( ((rank==0) && (r==0)) || ((rank==numprocs-1) && (r==rsize)) || c==0 || c==width-1 )
				sum = 0;
				
			else{
				for(x=-1; x<2; x++){
					for(y=-1; y<2; y++){
						sum = sum + submat[((r+row+x)*w)+(c+y)] * mask[x+1][y+1];
					}
				}
				if(sum < 0)
					sum = 0;
				else if(sum > 255)
					sum = 255;
				A[(r*w)+c] = sum;

			}

			
		} // end of column traversal
		
		
	} // end of row traversal


	time = MPI_Wtime() - time;
	//printf("Gathering \n");

	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	//printf("Gathered \n");
	
	Totalroutine[index] += time;



	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A); free(submat);

	return 1;
}

/*=======================================================================================*/
/*======================================= SOBEL =========================================*/
/*=======================================================================================*/

int edge_sobel(unsigned char *matrix, int w, int comp){

	
	double time;
	int row = 0;
	int i, j, x, y, sumx, sumy, shift;
	unsigned long r, c, rsize, rsize2, rowsize;

	// Limited to int by the MPI function call
	unsigned int disp[numprocs], disp2[numprocs], nsent[numprocs], nsent2[numprocs];
	unsigned long location = 0;

	int maskx[3][3];
	int masky[3][3];


	maskx[0][0] =-1; maskx[0][1] = 0; maskx[0][2] = 1;
	maskx[1][0] =-2; maskx[1][1] = 0; maskx[1][2] = 2;
	maskx[2][0] =-1; maskx[2][1] = 0; maskx[2][2] = 1;

	masky[0][0] =-1; masky[0][1] =-2; masky[0][2] =-1;
	masky[1][0] = 0; masky[1][1] = 0; masky[1][2] = 0;
	masky[2][0] = 1; masky[2][1] = 2; masky[2][2] = 1;


	rem = height % numprocs;
	quo = height / numprocs;
	start = quo * rank;

	// Load balancing
	if(rank < rem)
		rsize = quo + 1;
	else
		rsize = quo;


	if((rank == 0) || (rank == numprocs-1))
		rsize2 = rsize + 1; // rank 0 has lower row, and last rank has upper row;
	else
		rsize2 = rsize + 2; // middle ranks has upper and lower rows;


	unsigned char *A = (unsigned char *)malloc(rsize*w*comp*sizeof(unsigned char));
	unsigned char *submat = (unsigned char *)malloc(rsize2*w*comp*sizeof(unsigned char));


	// set the buffer size
	unsigned long bsize = rsize * w;
	unsigned long bsize2 = rsize2 * w;

	shift = 0;
	if(rank == 0){
		// Values for A[][]
		for(i=0;i<numprocs;i++){
			if(i<rem)
				rowsize = quo + 1;
			else
				rowsize = quo;
			nsent[i] = rowsize * w;
			disp[i] = i * quo * w;
			if(i<=rem){
				disp[i] = disp[i] + (i * w);
				if(i == rem)
					shift = (i * w);
			}
			else
				disp[i] = disp[i] + shift;
		}
		// Values for submat[][]
		for(i=0;i<numprocs;i++){
			if(i!=0)
				disp2[i] = disp[i] - w;
			else
				disp2[i] = disp[i];
			if(i == numprocs-1)
				nsent2[i] = nsent[i] + w;
			else{
				nsent2[i] = nsent[i] + (2 * w);
				if(i==0)
					nsent2[i] = nsent[i] + w;
			}
		}
		

	} // end of root task



	if(rank != 0)
		row = 1;

	//MPI_Barrier(MPI_COMM_WORLD);

	//printf("Scattering \n");
	
	MPI_Scatterv(matrix,nsent,disp,MPI_UNSIGNED_CHAR,A,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	MPI_Scatterv(matrix,nsent2,disp2,MPI_UNSIGNED_CHAR,submat,bsize2,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);	
	
	//printf("scattered \n");

	time = MPI_Wtime();


	for(r=0; r<rsize; r++){
		for(c=0; c<width; c++){

			int sum = 0;
			sumx = 0; sumy = 0;
			if( ((rank==0) && (r==0)) || ((rank==numprocs-1) && (r==rsize)) || c==0 || c==width-1 )
				sum = 0;
				
			else{
				for(x=-1; x<2; x++){
					for(y=-1; y<2; y++){
						sumx = sumx + submat[((r+row+x)*w)+(c+y)] * maskx[x+1][y+1];
					}
				}
				if(sumx < 0)
					sumx = 0;
				else if(sumx > 255)
					sumx = 255;


				for(x=-1; x<2; x++){
					for(y=-1; y<2; y++){
						sumy = sumy + submat[((r+row+x)*w)+(c+y)] * masky[x+1][y+1];
					}
				}
				if(sumy < 0)
					sumy = 0;
				else if(sumy > 255)
					sumy = 255;

				sum = abs(sumx) + abs(sumy);
				A[(r*w)+c] = sum;

			}

		} // end of column traversal
		
	} // end of row traversal

	time = MPI_Wtime() - time;

	
	MPI_Gather(A,bsize,MPI_UNSIGNED_CHAR,matrix,bsize,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);

	

	Totalroutine[9] += time;



	if(rank == 0){
	printf("Wall clock time: %f \n\n", time);
	//printf("Root access for raw image. \n");
		for(r=0;r<height;r++){ 
			for(i=0;i<w*comp;i++){ 
				//printf("width traversal \n");
				raw_image[location++] = matrix[r*w+i];
			}
		}
	}

	free(A); free(submat);


	return 1;
}



/**
 * write_jpeg_file Writes the raw image data stored in the raw_image buffer
 * to a jpeg image with default compression and smoothing options in the file
 * specified by *filename.
 *
 * \returns positive integer if successful, -1 otherwise
 * \param *filename char string specifying the file name to save to
 *
 */
int write_jpeg_file( char *filename )
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	
	/* this is a pointer to one row of image data */
	JSAMPROW row_pointer[1];
	FILE *outfile = fopen( filename, "wb" );
	
	if ( !outfile )
	{
		printf("Error opening output jpeg file %s! \n", filename );
		return -1;
	}
	cinfo.err = jpeg_std_error( &jerr );
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, outfile);

	/* Setting the parameters of the output file here */
	cinfo.image_width = width;	
	cinfo.image_height = height;
	cinfo.input_components = 1; //bytes_per_pixel;
	cinfo.in_color_space = cspace; //color_space;
    /* default compression parameters, we shouldn't be worried about these */
	jpeg_set_defaults( &cinfo );
	/* Now do the compression .. */
	jpeg_start_compress( &cinfo, TRUE );

	int component;
	if(icomp == 1)
		component = 1;
	else
		component = 3;

	/* like reading a file, this time write one row at a time */
	while( cinfo.next_scanline < cinfo.image_height )
	{
		row_pointer[0] = &raw_image[ cinfo.next_scanline * cinfo.image_width * component];
		jpeg_write_scanlines( &cinfo, row_pointer, 1 );
	}

	printf( "JPEG File Information after transformation: \n" );
	printf( "Image width and height: %d pixels and %d pixels.\n", cinfo.image_width, cinfo.image_height );
	printf( "Color components per pixel: %d.\n", cinfo.num_components );
	printf( "Color space: %d.\n\n", cinfo.jpeg_color_space );

	/* similar to read file, clean up after we're done compressing */
	jpeg_finish_compress( &cinfo );
	jpeg_destroy_compress( &cinfo );
	fclose( outfile );
	/* success code is 1! */
	return 1;
}

int main(int argc, char *argv[]){

//system("clear");
int rcount = 0;
int tcount = 0;
int i;

  //char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  

  //MPI_Get_processor_name(processor_name, &namelen);


  //printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

	char *infilename = "test.jpg", *outfilename = "test_out.jpg";

	infilename = argv[1];

	// Open a JPEG file
	if( read_jpeg_file( infilename ) > 0 ) 
	{
		// Create a new container file for the edited JPEG
		if(rank == 0){
			printf("Going to write the data to file. \n");
			if( write_jpeg_file( outfilename ) < 0 ) return -1;
		}
	}
	else return -1;

	if(rank == 0){
		printf("Time information: \n");
		for(i=0;i<rout;i++){
			printf("[%s] \nTotal exection time: %f \tAverage exec time: %f \n",
				routines[i],Totalroutine[i],Avgroutine[i]);
			if(Totalroutine[i] > 0){
				tcount++;
			}
			rcount += rcnt[i];
			Exectime += Totalroutine[i];
		}
		printf("\nOverall execution time: %f \n",Exectime);
		if(rcount > 0)
		AvgExec = Exectime / rcount;
		printf("Average execution time per any routine: %f \n",AvgExec);
		if(rcount > 0)
		AvgExec = Exectime / tcount;
		printf("Average execution time by routine type: %f \n\n",AvgExec);
	}

  MPI_Finalize();
return 0;
}


