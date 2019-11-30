#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<time.h>

typedef struct
{
    int * array;
    int index;
} Bucket;

void bucket_insert(Bucket * b, int x)
{
    b->array[b->index] = x;
    b->index = b->index +1;
}

int compare( const void * n1, const void * n2)
{
    return (*(int*)n1 - *(int*)n2);
}

// Checking for sortedness of a bucket
int check_bucket(Bucket * b)
{
    int i, last = b->array[0];
    for(i = 0; i < b->index; i++)
    {
	if(b->array[i] < last)
	    return -1;
	last = b->array[i];
    }
    return 0;
}


int main(int argc, char *argv[])
{
    int rank, size = 10, array_size = 100, max_num = 1000;    

    if(argc != 3){
        printf("\nERROR : Args reuired <array_size> <max_limit>");
        exit(1);
    }else{
        array_size = atio(argv[1]);
        max_num = atoi(argv[2]);
    }

    MPI_Init(&argc, &argv)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start, end;        
    int i, num_buckets, large_bucket_size, small_bucket_size;

    num_buckets = size;
    large_bucket_size = ceil(array_size/(float)size);
    small_bucket_size = large_bucket_size;
    int * array;

    if(rank == 0)
    {
        start = MPI_Wtime();
        array = (int *) malloc(sizeof(int) * array_size);
        srand(time(NULL));
        
        /* Initialize array with random numbers */
        for(i = 0; i < array_size; i++)
            array[i] = rand()%(max_num);

        printf(stdout,"No of buckets: %d\n Large bucket size:\t%d\n Small bucket size:\t%d\n Array size:\t\t%d\n ", 
            num_buckets,large_bucket_size,small_bucket_size,array_size);
    }

    //small buckets
    Bucket ** buckets = (Bucket **) malloc(sizeof(Bucket*)*num_buckets);
    
    for(i = 0; i < num_buckets; i++)
    {
        buckets[i] = (Bucket *) malloc(sizeof(Bucket));
        buckets[i]->array = (int*) malloc(sizeof(int)*small_bucket_size*2.0);
        buckets[i]->index = 0;
    }

    Bucket large_bucket;

    int * my_bucket_array = (int*) malloc(sizeof(int)*large_bucket_size*4.0);
    large_bucket.array = my_bucket_array;
    large_bucket.index = 0;
        
    /* Scatter array into large buckets */
    MPI_Scatter(array,large_bucket_size,MPI_INT,large_bucket.array,large_bucket_size,MPI_INT,0,MPI_COMM_WORLD);
    
    int dest;
    
    //divding no into samller buckets
    for(i = 0; i < large_bucket_size; i++)
    {
        dest = (large_bucket.array[i] * num_buckets)/max_num;
       
        if(dest == rank)
            bucket_insert(&large_bucket,large_bucket.array[i]);
        else
            bucket_insert(buckets[dest],large_bucket.array[i]);
    }

    MPI_Request * requests = (MPI_Request *) malloc(sizeof(MPI_Request) * size);
    
    for(i = 0; i < size; i++)
       if(i != rank)
           MPI_Isend(buckets[i]->array,small_bucket_size*2,MPI_INT,i,buckets[i]->index,MPI_COMM_WORLD,&requests[i]);
    
    MPI_Status status;
    int current = large_bucket.index;
    
    for(i = 0; i < size-1; i++)
    {
        MPI_Recv(&large_bucket.array[current],small_bucket_size*2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        current += status.MPI_TAG;
    }

    large_bucket.index = current;
  
    
    qsort(&large_bucket.array[0],current,sizeof(int),compare);
    
    // Get array sizes
    int * sizes = (int *) malloc(sizeof(int)*size);
    MPI_Gather(&current,1,MPI_INT,sizes,1,MPI_INT,0,MPI_COMM_WORLD);
    
    // Determine the displacement
    int * disp = (int *) malloc(sizeof(int)*size);
    if(rank == 0)
    {
        disp[0] = 0;
        for(i = 1; i < size+1; ++i)
            disp[i] = disp[i-1] + sizes[i-1];
    }

    /* Gathering results */
    MPI_Gatherv(large_bucket.array,current,MPI_INT,array,sizes,disp,MPI_INT,0,MPI_COMM_WORLD);
    
    if(rank == 0)
    {
        end = MPI_Wtime();
        printf("\narray size: %d\telapse time: %f\n",array_size,end-start);
        free(array);
    }

    MPI_Finalize();
    return 0;

}