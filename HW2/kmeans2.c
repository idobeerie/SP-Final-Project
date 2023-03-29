#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

typedef struct cord
{
    double value;
    struct cord *next;
}cord;

typedef struct vector
{
    struct vector *next;
    struct cord *cords;
    struct cord *prev_cords;
    int member_count;
    int choice;
}vector;



double distance(cord *cord1, cord *cord2){
    double d =0;
    cord *cptr1 = cord1;
    cord *cptr2 = cord2;
    while(cptr1 != NULL){
        d += pow((*cptr1).value - (*cptr2).value, 2);
        cptr1 = (*cptr1).next;
        cptr2 = (*cptr2).next;
    } 
    d = sqrt(d);
    return d;
}

void updateCord(cord *current, cord *newMember, int membersNum){
    double sum=0;
    cord *cptr1 = current;
    cord *cptr2 = newMember;
    while(cptr1 != NULL){
        sum = (*cptr1).value * membersNum + (*cptr2).value;
        (*cptr1).value = sum / (membersNum+1);
        cptr1 = (*cptr1).next;
        cptr2 = (*cptr2).next;
    }
}

void freeList(cord* head)
{
   cord* tmp;

   while (head != NULL)
    {
       tmp = head;
       head = head->next;
        free(tmp);
    }

}

void freeVector(vector* head)
{
   vector* temp;

   while (head != NULL)
    {
        freeList(head->cords);
        freeList(head->prev_cords);
        temp = head;
        head = head->next;
        free(temp);
    }

}


int main(int argc, char *argv[])
{
    vector *head_vec, *curr_vec, *head_d_vec, *curr_d_vec;
    cord *head_cord, *curr_cord, *head_d_cord, *curr_d_cord, *prev_cord;
    int iter = 200, k=0, count =0, i , min_i=0, iter_count=0, line_count=0, line_length=0, j=0, small_change=0, m=0;
    double epsi = 0.001, min_d = DBL_MAX, num=0.0, d=0.0;
    char line[1000];
    char *txt;
    FILE *fp;



    head_cord = malloc(sizeof(cord));
    curr_cord = head_cord;
    curr_cord->next = NULL;
    curr_cord->value =0;


    head_vec = malloc(sizeof(vector));
    curr_vec = head_vec;
    curr_vec->next = NULL;
    curr_vec->cords = NULL;
    curr_vec->prev_cords = NULL;
    curr_vec->member_count = 0;
    curr_vec->choice = 0;

    head_d_cord = malloc(sizeof(cord));
    curr_d_cord = head_d_cord;
    curr_d_cord->next = NULL;
    curr_d_cord->value = 0;

    head_d_vec = malloc(sizeof(vector));
    curr_d_vec = head_d_vec;
    curr_d_vec->next = NULL;
    curr_d_vec->cords = NULL;
    curr_d_vec->prev_cords = NULL;
    curr_d_vec->member_count = 0;
    curr_d_vec->choice = 0;



    for(i =0; i < k; i++){
        if(fp != NULL){
            fgets(line, sizeof(line), fp);
            txt = strtok(line, ",");
        }
        while(txt != NULL){
            if(txt != NULL)
            {
                num = atof(txt);
                curr_cord->value = num;
                curr_cord->next = malloc(sizeof(cord));
                curr_cord = curr_cord->next;
                curr_cord->next = NULL;
                curr_cord->value = 0;
                curr_cord->next = NULL;
            }
            line_length++;
            txt = strtok(NULL, ",");
        }
        curr_vec->member_count = 1;
        curr_vec->cords = head_cord;
        curr_vec->next = malloc(sizeof(vector));
        curr_vec->prev_cords = malloc(sizeof(cord));
        curr_vec = curr_vec->next;
        curr_vec->next = NULL;
        curr_vec->cords = NULL;
        curr_vec->prev_cords = NULL;
        curr_vec->member_count = 0;
        curr_vec->choice = 0;
        head_cord = malloc(sizeof(cord));
        curr_cord = head_cord;
        curr_cord->next = NULL;
        curr_cord->value = 0;
        curr_vec->cords = curr_cord;
        count++;
    }
    line_length = line_length / k;
    if(iter <= 1 || iter >= 1000){
        printf("Invalid maximum iteration!");
        exit(1);
    }
    rewind(fp);
    
    while(fgets(line, sizeof(line), fp)){
        line_count++;
        txt = strtok(line, ",");
        while(txt != NULL){
            if(txt != NULL){
                num = atof(txt);

                curr_d_cord->value = num;
                curr_d_cord->next = malloc(sizeof(struct cord));
                curr_d_cord = curr_d_cord->next;
                curr_d_cord->next = NULL;
                curr_d_cord->value = 0;
            }
            txt = strtok(NULL, ",");
        }
        curr_d_vec->cords = head_d_cord;
        curr_d_vec->next = malloc(sizeof(vector));
        curr_d_vec = curr_d_vec->next;
        curr_d_vec->next = NULL;
        curr_d_vec->cords = NULL;
        curr_d_vec->prev_cords = NULL;
        curr_d_vec->member_count = 0;
        curr_d_vec->choice = 0;
        head_d_cord = malloc(sizeof(cord));
        curr_d_cord = head_d_cord;
        curr_d_cord->next = NULL;
        curr_d_cord->value = 0;
        curr_d_vec->cords = curr_d_cord;
    }
    
    
    if(k >= line_count-1 || k <= 1){
        printf("Invalid number of clusters!");
        exit(1);
    }

    count =0;
    while(iter_count < iter){
        curr_d_vec = head_d_vec;
        curr_d_cord = curr_d_vec->cords;
        curr_vec = head_vec;
        curr_cord = head_vec->cords;
        while(count < line_count){
            min_d = DBL_MAX-1;
            for(i=0; i < k; i++){
                d = distance(curr_cord, curr_d_cord);
                if(d < min_d){
                    min_d = d;
                    min_i = i;
                }
                curr_vec = curr_vec->next;
                curr_cord = curr_vec->cords;
            }
            curr_d_vec->choice = min_i;            
            curr_vec = head_vec;
            curr_cord = head_vec->cords;
            curr_vec = head_vec;
            curr_cord = curr_vec->cords;
            curr_d_vec = curr_d_vec->next;
            curr_d_cord = curr_d_vec->cords;
            count++;
        }
        count = 0;
        curr_vec = head_vec;
        curr_cord = curr_vec->cords;
        prev_cord = curr_vec->prev_cords;
        for(i =0; i < k; i++){
            for(m =0; m < line_length; m++){
                prev_cord->value = curr_cord->value;
                if(iter_count == 0){
                    prev_cord->next = malloc(sizeof(cord));
                }
                prev_cord = prev_cord->next;
                prev_cord->value = 0;
                if(iter_count == 0){
                    prev_cord->next = NULL;
                }
                curr_cord = curr_cord->next;
            }
            curr_vec = curr_vec->next;
            curr_cord = curr_vec->cords;
            prev_cord = curr_vec->prev_cords;
        }
        curr_d_vec = head_d_vec;
        curr_vec = head_vec;
        for(i=0; i < k; i++){
            curr_cord = curr_vec->cords;
            while(curr_cord != NULL){
                curr_cord->value =0;
                curr_cord = curr_cord->next;
            }
            curr_vec->member_count = 0;
            curr_vec = curr_vec->next;
        }
        curr_d_vec = head_d_vec;
        while(count < line_count){
            curr_vec = head_vec;
            for(i=0; i < curr_d_vec->choice; i++){
                curr_vec = curr_vec->next;            
            }
            updateCord(curr_vec->cords, curr_d_vec->cords, curr_vec->member_count);
            curr_vec->member_count+=1;
            curr_d_vec = curr_d_vec->next;
            count++;       
        }
        count = 0;
        curr_vec = head_vec;
        if(iter_count > 0){
            small_change = 1;
            for(i =0; i < k; i++){
                if(distance(curr_vec->cords, curr_vec->prev_cords) >= epsi){
                    small_change = 0;
                }
                curr_vec = curr_vec->next;
            }
        }
        if(small_change == 1){
            break;
        }
        iter_count++;
    }

    curr_vec = head_vec;
    curr_cord = (*curr_vec).cords;
    for(i=0; i < k; i++){
        curr_cord = (*curr_vec).cords;
        for(j=0; j < line_length; j++){
            printf("%.4f, ", (*curr_cord).value);
            curr_cord = (*curr_cord).next;
        }
        printf("\n");
        curr_vec = (*curr_vec).next;
    }
    fclose(fp);
    freeVector(head_vec);
    freeVector(head_d_vec);
    return 0;   
}