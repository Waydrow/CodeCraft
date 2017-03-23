#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

typedef struct Chrom {                         // �ṹ�����ͣ�Ϊ����Ⱦɫ��Ľṹ��
    short int bit[6];//һ��6bit����Ⱦɫ����б��룬����1λΪ����λ��ȡֵ��Χ-64~+64
    int fit ;//��Ӧֵ
    double rfit;//��Ե�fitֵ������ռ�İٷֱ�
    double cfit;//���۸���
} chrom;
//���彫���õ��ļ���������
void *evpop (chrom popcurrent[4]);//������Ⱥ�ĳ�ʼ��
int x (chrom popcurrent);
int y (int x);
void *pickchroms (chrom popcurrent[4]);//ѡ�����
void *pickchroms_new (chrom popcurrent[4]); // ���ڸ��ʷֲ�
void *crossover (chrom popnext[4]);//�������
void *mutation (chrom popnext[4]);//ͻ��
double r8_uniform_ab ( double a, double b, int &seed );//����a~b֮����ȷֲ�������
chrom popcurrent [4];                        // ��ʼ��Ⱥ��ģΪ��
chrom popnext [4];                           // ���º���Ⱥ��ģ��Ϊ��
int main () {                                  // ��������
    int num ;                                    // ����������
    int i ,j, l,Max ,k;
    Max=0;                                      // �������ֵ

    printf("\nWelcome to the Genetic Algorithm��\n");  //
    printf("The Algorithm is based on the function y = -x^2 + 5 to find the maximum value of the function.\n");

enter:
    printf ("\nPlease enter the no. of iterations\n��������Ҫ�趨�ĵ����� : ");
    scanf("%d" ,&num);                           // ����������������͸����� num��

    if(num <1)
        goto enter ;                                 // �ж�����ĵ��������Ƿ�Ϊ�����㣬�ǵĻ��������룻
    //��ͬ����������ܽ����ͬ�������ǵ������õĵ�����������ʱ��Ⱦɫ��Ļ����͹��������ֲ�����
    srand(time(0));
    evpop(popcurrent );    // ���������ʼ��Ⱥ��
    //�Ƿ���Ҫָ��x��ȡֵ��Χ�أ�6bit����ʾ���֣���һλΪ����λ��5bit��ʾ���ִ�С�����ԣ�ȡֵ��ΧΪ-32~+31
    Max = popcurrent[0].fit;//��Maxֵ���г�ʼ��

    for(i =0; i< num; i ++) {                      // ��ʼ������

        printf("\ni = %d\n" ,i);                 // �����ǰ����������

        for(j =0; j<4; j++) {
            popnext[j ]=popcurrent[ j];           // ������Ⱥ��
        }

        pickchroms(popnext );                    // ��ѡ������壻
        crossover(popnext );                     // ����õ��¸��壻
        mutation(popnext );                      // ����õ��¸��壻

        for(j =0; j<4; j++) {
            popcurrent[j ]=popnext[ j];              // ��Ⱥ���棻
        }

    }  // �ȴ�������ֹ��
//�����������������Ҫע��ȡ�ϴ�ĵ�������
    for(l =0; l<3; l++) {
        if(popcurrent [l]. fit > Max ) {
            Max=popcurrent [l]. fit;
            k=x(popcurrent [l]);//��ʱ��value��Ϊ�����xֵ
        }

    }
    printf("\n ��x���� %dʱ�������õ����ֵΪ�� %d ",k ,Max);
    printf("\nPress any key to end ! " );

    _flushall();                                 // ������л�������
    getche();                                   // �ӿ���̨ȡ�ַ������Իس�Ϊ������
    return 0;

}



void *evpop (chrom popcurrent[4]) { // ������������ɳ�ʼ��Ⱥ��
    int i ,j, value1;
    int random ;
    double sum=0;

    for(j =0; j<4; j++) {                         // ����Ⱥ�еĵ�1��Ⱦɫ�嵽��4��Ⱦɫ��
        for(i =0; i<6; i++) {                    // ��Ⱦɫ��ĵ�1������λ����6������λ
            random=rand ();                     // ����һ�����ֵ
            random=(random %2);                 // �������0����1
            popcurrent[j ].bit[ i]=random ;       // �������Ⱦɫ����ÿһ������λ��ֵ����
        }

        value1=x (popcurrent[ j]);                // �������ƻ���Ϊʮ���ƣ��õ�һ������ֵ��
        popcurrent[j ].fit= y(value1); // ����Ⱦɫ�����Ӧ��ֵ
        sum = sum + popcurrent[j ].fit;
        printf("\n popcurrent[%d]=%d%d%d%d%d%d  value=%d  fitness = %d",j, popcurrent[j ].bit[5], popcurrent[j ].bit[4], popcurrent[j ].bit[3], popcurrent[j ].bit[2], popcurrent[j ].bit[1], popcurrent[j ].bit[0], value1,popcurrent [j]. fit);
        // �������Ⱦɫ��ı��������
    }
    //������Ӧֵ�ðٷֱȣ��ò������������̶�ѡ��ʱ��Ҫ�õ���
    for (j = 0; j < 4; j++) {
        popcurrent[j].rfit = popcurrent[j].fit/sum;
        popcurrent[j].cfit = 0;//�����ʼ��Ϊ0
    }
    return(0);
}


int x (chrom popcurrent) { // �������������ƻ���Ϊʮ���ƣ�
    //�˴���Ⱦɫ�峤��Ϊ�����и���ʾ����λ

    int z ;
    z=(popcurrent .bit[0]*1)+( popcurrent.bit [1]*2)+(popcurrent. bit[2]*4)+(popcurrent .bit[3]*8)+( popcurrent.bit [4]*16);

    if(popcurrent .bit[5]==1) { // ���ǵ����ţ�
        z=z *(-1);
    }

    return(z );
}
//��Ҫ���ܹ����ⲿֱ�Ӵ��亯������ǿ³����
int y (int x) { // ��������������Ӧ�ȣ�
    int y ;
    y=-(x *x)+5;                                // Ŀ�꺯����y= - ( x^ 2 ) +5��
    return(y );
}
//�������̶�ѡ�񷽷������л����͵�ѡ��
void *pickchroms_new (chrom popnext[4]) { //�������
    int men;
    int i;
    int j;
    double p;
    double sum=0.0;
    //find the total fitness of the population
    for (men = 0; men < 4; men++ ) {
        sum = sum + popnext[men].fit;
    }
    //calculate the relative fitness of each member
    for (men = 0; men < 4; men++ ) {
        popnext[men].rfit = popnext[men].fit / sum;
    }
    //calculate the cumulative fitness,��������۸���
    popcurrent[0].cfit = popcurrent[0].rfit;
    for ( men = 1; men < 4; men++) {
        popnext[men].cfit = popnext[men-1].cfit + popnext[men].rfit;
    }

    for ( i = 0; i < 4; i++ ) {
        //����0~1֮��������
        //p = r8_uniform_ab ( 0, 1, seed );//ͨ����������0~1֮����ȷֲ�������
        p =rand()%10;//
        p = p/10;
        if ( p < popnext[0].cfit ) {
            popcurrent[i] = popnext[0];
        } else {
            for ( j = 0; j < 4; j++ ) {
                if ( popnext[j].cfit <= p && p < popnext[j+1].cfit ) {
                    popcurrent[i] = popcurrent[j+1];
                }
            }
        }
    }
    //  Overwrite the old population with the new one.
    //
    for ( i = 0; i < 4; i++ ) {
        popnext[i] = popcurrent[i];
    }
    return(0);
}
void *pickchroms (chrom popnext[4]) {        // ������ѡ����壻
    int i ,j;
    chrom temp ;                                // �м����
    //��˴˴���Ƶ��Ǹ����壬���Բ�����
    for(i =0; i<3; i++) {                        // ���ݸ�����Ӧ�������򣻣�ð�ݷ���
        for(j =0; j<3-i; j++) {
            if(popnext [j+1]. fit>popnext [j]. fit) {
                temp=popnext [j+1];
                popnext[j +1]=popnext[ j];
                popnext[j ]=temp;

            }
        }
    }
    for(i =0; i<4; i++) {
        printf("\nSorting:popnext[%d] fitness=%d" ,i, popnext[i ].fit);
        printf("\n" );
    }
    _flushall();/* ������л����� */
    return(0);
}
double r8_uniform_ab( double a, double b, int &seed ) {
    {
        int i4_huge = 2147483647;
        int k;
        double value;

        if ( seed == 0 ) {
            std::cerr << "\n";
            std::cerr << "R8_UNIFORM_AB - Fatal error!\n";
            std::cerr << "  Input value of SEED = 0.\n";
            exit ( 1 );
        }

        k = seed / 127773;

        seed = 16807 * ( seed - k * 127773 ) - k * 2836;

        if ( seed < 0 ) {
            seed = seed + i4_huge;
        }

        value = ( double ) ( seed ) * 4.656612875E-10;

        value = a + ( b - a ) * value;

        return value;
    }
}
void *crossover (chrom popnext[4]) {            // ���������������

    int random ;
    int i ;
    //srand(time(0));
    random=rand ();                             // �����������㣻
    random=((random %5)+1);                     // ����������0��5֮�䣻
    for(i =0; i< random; i ++) {
        popnext[2].bit [i]= popnext[0].bit [i];   // child 1 cross over
        popnext[3].bit [i]= popnext[1].bit [i];   // child 2 cross over
    }

    for(i =random; i<6; i ++) {                   // crossing the bits beyond the cross point index
        popnext[2].bit [i]= popnext[1].bit [i];    // child 1 cross over
        popnext[3].bit [i]= popnext[0].bit [i];    // chlid 2 cross over
    }

    for(i =0; i<4; i++) {
        popnext[i ].fit= y(x (popnext[ i]));        // Ϊ�¸��������Ӧ��ֵ��
    }

    for(i =0; i<4; i++) {
        printf("\nCrossOver popnext[%d]=%d%d%d%d%d%d    value=%d    fitness = %d",i, popnext[i ].bit[5], popnext[i ].bit[4], popnext[i ].bit[3], popnext[i ].bit[2], popnext[i ].bit[1], popnext[i ].bit[0], x(popnext [i]), popnext[i ].fit);
        // ����¸��壻
    }
    return(0);
}

void *mutation (chrom popnext[4]) {             // ���������������

    int random ;
    int row ,col, value;
    //srand(time(0));
    random=rand ()%50;  // ���������֮�������
    //�������ҲҪ���һ���ĸ��������У�һ������Ϊ0��0.5֮��
    //
    if(random ==25) {                            // random==25�ĸ���ֻ��2%����������Ϊ����������С���ʽ��б��죡��
        col=rand ()%6;                            // �������Ҫ����Ļ���λ�ţ�
        row=rand ()%4;                            // �������Ҫ�����Ⱦɫ��ţ�

        if(popnext [row]. bit[col ]==0) {           // 1��Ϊ��
            popnext[row ].bit[ col]=1 ;
        } else if (popnext[ row].bit [col]==1) {    // 0��Ϊ��
            popnext[row ].bit[ col]=0;
        }
        popnext[row ].fit= y(x (popnext[ row]));     // �����������Ӧ��ֵ��
        value=x (popnext[ row]);
        printf("\nMutation occured in popnext[%d] bit[%d]:=%d%d%d%d%d%d    value=%d   fitness=%d", row,col ,popnext[ row].bit [5],popnext[ row].bit [4],popnext[ row].bit [3],popnext[ row].bit [2],popnext[ row].bit [1],popnext[ row].bit [0],value, popnext[row ].fit);

        // ����������¸��壻
    }

    return(0);
}
