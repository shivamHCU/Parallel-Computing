#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<mpi.h>
#include<math.h>

typedef struct
{
    int intervals,rank,size;
    double delta,region;
    double min,max;
} Params;

double f(double x)
{
    return sqrt(1-x*x); 
}

// f(x) using the rectangle rule
double rectangle_rule(double start, double end, double delta)
{
    double area = 0.0;
    double x;
    for(x = start; x <= end; x+= delta)
        area += f(x)*delta;
    return area;
}

void master(Params * p)
{   
    double last = 0.0;
    double area = p->error+1,local;
    double start, end;

    while(fabs(last - area) > p->error)
    {
        last = area;
        start = MPI_Wtime();
        MPI_Bcast(&p->intervals,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Reduce(&local,&area,1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        end = MPI_Wtime();
        p->delta = (p->max-p->min)/p->intervals;
        fprintf(stderr,"Intervals = %8d\tDelta = %f\tResult = %f\tError=%f\tTime = %f\n",p->intervals,p->delta,area,fabs(last-area),end-start);
        p->intervals = p->intervals*2;
    } 

    p->intervals = -1;
    MPI_Bcast(&p->intervals,1,MPI_INT,0,MPI_COMM_WORLD);
    return;
}

void slave(Params *p)
{
    double start, end, area;
    while(1)
    {
        MPI_Bcast(&p->intervals,1,MPI_INT,0,MPI_COMM_WORLD);
        if(p->intervals == -1)
            return;
        p->region = (p->max-p->min)/(double)p->size;
        start = p->min + p->region * (p->rank-1);
        end = start + p->region;
        p->delta = (p->max-p->min)/p->intervals;
        area = rectangle_rule(start,end,p->delta);
        MPI_Reduce(&area,&area,1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
}



int main(int argc, char *argv[])
{
    int rank, size;
    float min, max, intervals;

    if(argc != 4){
        printf("\nERROR : Args reuired <start> <end> <intervals>");
        exit(1);
    }else{
        min = atof(argv[1]);
        max = atof(argv[2]);
        intervals = atof(argv[3]);
    }

    if(MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        fprintf(stderr, "Unable to initialize MPI!\n");
        return -1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    Params p;
    p.min = min;
    p.max = max;
    p.error = 0.00001;
    p.intervals = intervals;
    p.rank = rank;
    p.size = size-1;
    
    if(rank == 0)
        master(&p);
    else
        slave(&p);

    MPI_Finalize();
    return 0;
}