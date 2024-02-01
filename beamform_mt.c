// 3D Ultrasound beamforming baseline code for EECS 570 
// Created by: Richard Sampson, Amlan Nayak, Thomas F. Wenisch
// Revision 1.0 - 11/15/16

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdint.h>
#include <pthread.h>

struct threadStruct {
	int threadID;
	int subsetLength;
	int sls_t;
	int sls_p;
	int pts_r;
	float tx_x;
	float tx_y;
	float tx_z;
	float *point_x; 
	float *point_y;
	float *point_z; 
	float* image;
	int offset;
	int data_len;
	float rx_z;
	float idx_const;
	int filter_delay;

};

float *dist_tx; // Transmit distance (ie first leg only)
float *image;
float *rx_x; // Receive transducer x position
float *rx_y; // Receive transducer y position
float *rx_data; // Pointer to pre-processed receive channel data
float *point_x; // Point x position
float *point_y; // Point y position
float *point_z; // Point z position
pthread_mutex_t mutexImage;


void *firstLoop(void *arg){
	struct threadStruct *myArg;
	myArg = (struct threadStruct *)arg;

	int i, start, end;
	float x_comp, y_comp, z_comp;
	//float *point_x = myArg->point_x;
	//float *point_y = myArg->point_y;
	//float *point_z = myArg->point_z;
	int ID = myArg->threadID;

	start = (myArg->threadID) * (myArg->subsetLength);
	end = start + (myArg->subsetLength); 

	//printf("Thread: %d, Start: %d, End: %d\n", ID, start, end);

	for(i = start; i < end; i++){
		x_comp = (myArg->tx_x) - point_x[i];
		x_comp = x_comp * x_comp;
		y_comp = (myArg->tx_y) - point_y[i];
		y_comp = y_comp * y_comp;
		z_comp = (myArg->tx_z) - point_z[i];
		z_comp = z_comp * z_comp;

        //if(i > 399360) printf("ARRAY ACCESS OUT OF BOUNDS#1!");

		dist_tx[i] = (float)sqrt(x_comp + y_comp + z_comp);

	}
    //printf("ACCESSING dist_tx out of bounds! \n");
    //dist_tx[400,000];
	pthread_exit(NULL);
}

void *secondLoop(void *arg){
	struct threadStruct *myArg;
	myArg = (struct threadStruct *)arg;

	int it_rx, index, start, end, i;
	float x_comp, y_comp, z_comp,dist;
	//float *point_x = myArg->point_x;
	//float *point_y = myArg->point_y;
	//float *point_z = myArg->point_z;
	int ID = myArg->threadID;
	int offset = 0;
	float *image_pos = image;
	int sls_t = myArg->sls_t;
	int sls_p = myArg->sls_p;
	int pts_r = myArg->pts_r;
	int data_len = myArg->data_len;
	float rx_z = myArg->rx_z;
	float idx_const = myArg->idx_const;
	int filter_delay = myArg->filter_delay;

	start = (myArg->threadID) * (myArg->subsetLength);
	end = start + (myArg->subsetLength);

	for (it_rx = start; it_rx < end; it_rx++) {

        // image_pos can't be global, or threads will increment each others image_pos
        //image_pos needs to be local and array must be returned to main
        //and then in main add up all the returned image_pos arrays tgth
		image_pos = image; // Reset image pointer back to beginning
		//point = 0; // Reset 
		//printf("it_rx = %d\n", it_rx);

		for(i = 0; i < (sls_t * sls_p * pts_r); i++) {
			
			x_comp = rx_x[it_rx] - point_x[i];
			x_comp = x_comp * x_comp;
			y_comp = rx_y[it_rx] - point_y[i];
			y_comp = y_comp * y_comp;
			z_comp = rx_z - point_z[i];
			z_comp = z_comp * z_comp;

			dist = dist_tx[i] + (float)sqrt(x_comp + y_comp + z_comp);
			index = (int)(dist/idx_const + filter_delay + 0.5);

			pthread_mutex_lock(&mutexImage);
			*image_pos++ += rx_data[index+offset];
			pthread_mutex_unlock(&mutexImage);

			//printf("i = %d\n", i);
		}

		offset += data_len;
	}
    //point_x = NULL;
    //point_y = NULL;
    //point_z = NULL;
    //image = NULL;
	pthread_exit(NULL);
}

int main (int argc, char **argv) {

	int size = atoi(argv[2]);

	/* Variables for transducer geometry */
	int trans_x = 32; // Transducers in x dim
	int trans_y = 32; // Transducers in y dim
	
	// float *rx_x; // Receive transducer x position
	// float *rx_y; // Receive transducer y position
	float rx_z = 0; // Receive transducer z position

	int data_len = 12308; // Number for pre-processed data values per channel
	int offset = 0; // Offset into rx_data
	// float *rx_data; // Pointer to pre-processed receive channel data

	float tx_x = 0; // Transmit transducer x position
	float tx_y = 0; // Transmit transducer y position
	float tx_z = -0.001; // Transmit transducer z position

	/* Variables for image space points */
	int point; // Index into image space

//	float *point_x; // Point x position
//	float *point_y; // Point y position
//	float *point_z; // Point z position

	int pts_r = 1560; // Radial points along scanline
	int sls_t = size; // Number of scanlines in theta
	int sls_p = size; // Number of scanlines in phi

	//float *image_pos; // Pointer to current position in image
	//float *image;  // Pointer to full image (accumulated so far)

	/* Iterators */
	int it_rx; // Iterator for recieve transducer
	int it_r; // Iterator for r
	int it_t; // Iterator for theta
	int it_p; // Iterator for phi

	/* Variables for distance calculation and index conversion */
	float x_comp; // Itermediate value for dist calc
	float y_comp; // Itermediate value for dist calc
	float z_comp; // Itermediate value for dist calc

	//float *dist_tx; // Transmit distance (ie first leg only)
	float dist; // Full distance
	const float idx_const = 0.000009625; // Speed of sound and sampling rate, converts dist to index
	const int filter_delay = 140; // Constant added to index to account filter delay (off by 1 from MATLAB)
	int index; // Index into transducer data

        FILE* input;
        FILE* output;

	/* Allocate space for data */
	rx_x = (float*) malloc(trans_x * trans_y * sizeof(float));
	if (rx_x == NULL) fprintf(stderr, "Bad malloc on rx_x\n");
	rx_y = (float*) malloc(trans_x * trans_y * sizeof(float));
	if (rx_y == NULL) fprintf(stderr, "Bad malloc on rx_y\n");
	rx_data = (float*) malloc(data_len * trans_x * trans_y * sizeof(float));
	if (rx_data == NULL) fprintf(stderr, "Bad malloc on rx_data\n");

	point_x = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (point_x == NULL) fprintf(stderr, "Bad malloc on point_x\n");
	point_y = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (point_y == NULL) fprintf(stderr, "Bad malloc on point_y\n");
	point_z = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (point_z == NULL) fprintf(stderr, "Bad malloc on point_z\n");

	dist_tx = (float*) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (dist_tx == NULL) fprintf(stderr, "Bad malloc on dist_tx\n");

	image = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (image == NULL) fprintf(stderr, "Bad malloc on image\n");
	memset(image, 0, pts_r * sls_t * sls_p * sizeof(float));

	/* validate command line parameter */
	if (argc < 2 || !(strcmp(argv[2],"16") || strcmp(argv[2],"32") || strcmp(argv[2],"64"))) {
	  printf("Usage: %s {16|32|64}\n",argv[0]);
	  fflush(stdout);
	  exit(-1);
	}

    int numThreads = atoi(argv[1]);
    printf("NUMTHREADS = %s \n", argv[1]);

	char buff[128];
        #ifdef __MIC__
	  sprintf(buff, "/beamforming_input_%s.bin", argv[2]);
        #else // !__MIC__
	  sprintf(buff, "/cad2/ece1755s/assignment1_data/beamforming_input_%s.bin", argv[2]);
        #endif

        input = fopen(buff,"rb");
	if (!input) {
	  printf("Unable to open input file %s.\n", buff);
	  fflush(stdout);
	  exit(-1);
	}	

	/* Load data from binary */
	fread(rx_x, sizeof(float), trans_x * trans_y, input); 
	fread(rx_y, sizeof(float), trans_x * trans_y, input); 

	fread(point_x, sizeof(float), pts_r * sls_t * sls_p, input); 
	fread(point_y, sizeof(float), pts_r * sls_t * sls_p, input); 
	fread(point_z, sizeof(float), pts_r * sls_t * sls_p, input); 

	fread(rx_data, sizeof(float), data_len * trans_x * trans_y, input); 
        fclose(input);

	printf("Beginning computation\n");
	fflush(stdout);

	/* Create threads and thread structure arguments */

	struct threadStruct threadStructArray[numThreads];
	pthread_t threads[numThreads];
	struct threadStruct threadStructArray2[numThreads];
	pthread_t threads2[numThreads];
	pthread_mutex_init(&mutexImage,NULL);


	/* get start timestamp */
 	struct timeval tv;
    	gettimeofday(&tv,NULL);
    	uint64_t start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
 
	/* --------------------------- COMPUTATION ------------------------------ */

	// Initialize structs and create threads
	for (int x = 0; x < numThreads; x++){
		threadStructArray[x].threadID = x;
		threadStructArray[x].subsetLength = (sls_t * sls_p * pts_r)/numThreads;
		threadStructArray[x].sls_t = sls_t;
		threadStructArray[x].sls_p = sls_p;
		threadStructArray[x].pts_r = pts_r;
		threadStructArray[x].tx_x = tx_x;
		threadStructArray[x].tx_y = tx_y;
		threadStructArray[x].tx_z = tx_z;
		threadStructArray[x].point_x = point_x;
		threadStructArray[x].point_y = point_y;
		threadStructArray[x].point_z = point_z;

		if (pthread_create(&threads[x], NULL, firstLoop, (void*)&threadStructArray[x]) != 0) {
			fprintf(stderr, "Error creating thread %d\n", x);
			exit(-1);
		}
	}


	/* First compute transmit distance */
	// point = 0;
	// for (it_t = 0; it_t < sls_t; it_t++) {

	// 	for (it_p = 0; it_p < sls_p; it_p++) {
	// 		for (it_r = 0; it_r < pts_r; it_r++) {

	// 			x_comp = tx_x - point_x[point];
	// 			x_comp = x_comp * x_comp;
	// 			y_comp = tx_y - point_y[point];
	// 			y_comp = y_comp * y_comp;
	// 			z_comp = tx_z - point_z[point];
	// 			z_comp = z_comp * z_comp;

	// 			dist_tx[point++] = (floavalgrind --tool=memcheck ./your_programpoint];
	// 			y_comp = y_comp * y_comp;
	// 			z_comp = tx_z - point_z[point];
	// 			z_comp = z_comp * z_comp;

	// 			dist_tx[point++] = (float)sqrt(x_comp + y_comp + z_comp);		
	// }

	// printf("point = %d\n", point);
	// printf("sls_t = %d \n", sls_t);
	// printf("sls_p = %d\n", sls_p);
	// printf("pts_r = %d\n", pts_r);
	// printf("////////////////////");

	// for(int a = 0; a < 20; a++) {
	// 	printf("MT: distx[%d] = %f\n", a, dist_tx[a]);
	// }

	// Main thread waits for threads to finish their first loops
    for (int y = 0; y < numThreads; y++) {
        if (pthread_join(threads[y], NULL) != 0) {
            fprintf(stderr, "Error joining thread %d\n", y);
            exit(-1);
        }
    }

	// Create new threads for second loop
	// Initialize structs and create threads for second loop
	for (int x = 0; x < numThreads; x++){
		threadStructArray2[x].threadID = x;
		threadStructArray2[x].subsetLength = (trans_x * trans_y)/numThreads;
		threadStructArray2[x].sls_t = sls_t;
		threadStructArray2[x].sls_p = sls_p;
		threadStructArray2[x].pts_r = pts_r;
		threadStructArray2[x].tx_x = tx_x;
		threadStructArray2[x].tx_y = tx_y;
		threadStructArray2[x].tx_z = tx_z;
		threadStructArray2[x].point_x = point_x;
		threadStructArray2[x].point_y = point_y;
		threadStructArray2[x].point_z = point_z;
		threadStructArray2[x].image = image;
		threadStructArray2[x].data_len = data_len;
		threadStructArray2[x].rx_z = rx_z;
		threadStructArray2[x].idx_const = idx_const;
		threadStructArray2[x].filter_delay = filter_delay;
        threadStructArray2[x].offset = offset;



		if (pthread_create(&threads2[x], NULL, secondLoop, (void*)&threadStructArray2[x]) != 0) {
            fprintf(stderr, "Error creating thread %d\n", x);
            exit(-1);
        }
	}

	// Main thread waits for threads to finish their second loops
	for (int y = 0; y < numThreads; y++) {
        if (pthread_join(threads2[y], NULL) != 0) {
            fprintf(stderr, "Error joining thread %d\n", y);
            exit(-1);
        }
    }

	/* Now compute reflected distance, find index values, add to image */
//	 for (it_rx = 0; it_rx < trans_x * trans_y; it_rx++) {
//
//	 	image_pos = image; // Reset image pointer back to beginning
//	 	point = 0; // Reset
//
//	 	//Iterate over entire image space
//	 	for (it_t = 0; it_t < sls_t; it_t++) {
//	 		for (it_p = 0; it_p < sls_p; it_p++) {
//	 			for (it_r = 0; it_r < pts_r; it_r++) {
//
//	 				x_comp = rx_x[it_rx] - point_x[point];
//	 				x_comp = x_comp * x_comp;
//	 				y_comp = rx_y[it_rx] - point_y[point];
//	 				y_comp = y_comp * y_comp;
//	 				z_comp = rx_z - point_z[point];
//	 				z_comp = z_comp * z_comp;
//
//	 				dist = dist_tx[point++] + (float)sqrt(x_comp + y_comp + z_comp);
//	 				index = (int)(dist/idx_const + filter_delay + 0.5);
//	 				*image_pos++ += rx_data[index+offset];
//	 			}
//	 		}
//	 	}
//	 	offset += data_len;
//	}

	// make the three innermost loops into one loop
	// for (it_rx = 0; it_rx < trans_x * trans_y; it_rx++) {

	// 	image_pos = image; // Reset image pointer back to beginning
	// 	point = 0; // Reset 

	// 	for(int b = 0; b < (sls_t * sls_p * pts_r); b++) {
			
	// 		x_comp = rx_x[it_rx] - point_x[point];
	// 		x_comp = x_comp * x_comp;
	// 		y_comp = rx_y[it_rx] - point_y[point];
	// 		y_comp = y_comp * y_comp;
	// 		z_comp = rx_z - point_z[point];
	// 		z_comp = z_comp * z_comp;

	// 		dist = dist_tx[point++] + (float)sqrt(x_comp + y_comp + z_comp);
	// 		index = (int)(dist/idx_const + filter_delay + 0.5);
	// 		*image_pos++ += rx_data[index+offset];
	// 	}

	// 	offset += data_len;
	// }
	/* --------------------------------------------------------------------- */

	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	uint64_t end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	uint64_t elapsed = end - start;

	printf("@@@ Elapsed time (usec): %lu\n", elapsed);
	printf("Processing complete.  Preparing output.\n");
	fflush(stdout);

	/* Write result to file */
	char* out_filename;
        #ifdef __MIC__
	  out_filename = "/home/micuser/beamforming_output.bin";
        #else // !__MIC__
	  out_filename = "beamforming_output.bin";
        #endif
        output = fopen(out_filename,"wb");
	fwrite(image, sizeof(float), pts_r * sls_t * sls_p, output); 
	fclose(output);

	printf("Output complete.\n");
	fflush(stdout);

	/* Cleanup */


    free(rx_x);
	free(rx_y);
	free(rx_data);
	printf("FREE DONE\n");
	free(point_x);
    //printf("FREE DONE\n");
    //printf("point_y: %f\n", *(point_y+2));
    printf("XXXXXXXXXXXXXXX\n");
    free(point_y);
    printf("YYYYYYYYYYYYYYYY\n");
	free(point_z);
	printf("ZZZZZZZZZZZZZZZZZ\n");
	free(dist_tx);
    printf("DISTTTTTTTTTTTTTTT\n");
	free(image);
    printf("IMGGGGGGGGGGGGGGGG\n");
	pthread_mutex_destroy(&mutexImage);
    printf("MUTEEEEEEEEEEEEEEE\n");

	return 0;
}
