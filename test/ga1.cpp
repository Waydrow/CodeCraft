//�Ŵ��㷨 GA
#include<iostream>
#include <cstdlib>
#include<bitset>
using namespace std;
const int L=5; //�������ĳ���
int f(int x) { //������躯��f(x)
    int result;
    result=x*x*x-60*x*x+900*x+100;
    return result;
}
int main(int argc,char *argv[]) {
    int a(0),b(32); //����x�Ķ�����Χ
    const int pop_size=8; //������Ⱥ��С
// int L; //ָ������ĳ���
    const int NG=20; //ָ����Ⱥ���ķ�ֳ�Ĵ���
    int t=0; //��ǰ��ֳ�Ĵ���
    int p[pop_size]; //������Ⱥ
    int q[pop_size]; //���己ֳ��Ⱥ ����Ⱥ����һ��
    srand(6553); //������������ɵ�����
    double sum; //��ֵ�ܺ�
    double avl_sum; //�ʶ�ƽ��ֵ
    double p_probability[pop_size]; //��ֵ����
    double pp[pop_size];
    double pro; //����������ɵĸ���
    float pc=0.90; //���彻��ĸ���
    float pm=0.05; //�������ĸ���
    cout<<"��ʼ����Ⱥ ";
    for(int i=0; i<pop_size; i++) { //���ɳ�ʼ�ĵ�0����Ⱥ
        p[i]=rand()%31;
        cout<<p[i]<<" ";
    }
    cout<<endl;
    cout<<endl;
    void Xover(int &,int &); //�������溯��
//��ֹͣ׼������ ����ֳ����û�������� ,������ֳ
    while(t<=NG) {
        cout<<"��ֳ�Ĵ�����t="<<t<<endl;
        sum=0.0;
        for(int i=0; i<pop_size; i++) {
            q[i]=p[i];
            cout<<q[i]<<" ";
        }
        cout<<endl;
        for(int i=0; i<pop_size; i++) //����sum
            sum +=f(p[i]);
        avl_sum=sum/pop_size;
        cout<<"sum="<<sum<<endl;
        cout<<"�ʶ�ƽ��ֵ="<<avl_sum<<endl;
        for(int i=0; i<pop_size; i++) { //������ֵ����
            p_probability[i]=f(p[i])/sum;
            if(i==0) {
                pp[i]=p_probability[i];
                cout<<"pp"<<i<<"="<<pp[i]<<endl;
            } else {
                pp[i]=p_probability[i]+pp[i-1];
                cout<<"pp"<<i<<"="<<pp[i]<<endl;
            }
//cout<<"p_probability"<<i<<"="<<p_probability[i]<<endl;
        }
        //ѡ��˫��
        for(int i=0; i<pop_size; i++) {
            pro=rand()%1000/1000.0;
            if(pro>=pp[0]&&pro<pp[1])
                p[i]=q[0];
            else if(pro>=pp[1]&&pro<pp[2])
                p[i]=q[1];
            else if(pro>=pp[2]&&pro<pp[3])
                p[i]=q[2];
            else if(pro>=pp[3]&&pro<pp[4])
                p[i]=q[3];
            else if(pro>=pp[4]&&pro<pp[5])
                p[i]=q[4];
            else
                p[i]=q[5];
        }
//�ӽ�����
        int r=0;
        int z=0;
        for(int j=0; j<pop_size; j++) {
            pro=rand()%1000/1000.0;
            if(pro<pc) {
                ++z;
                if(z%2==0)
                    Xover(p[r],p[j]);
                else
                    r=j;
            }
        }
//��������
        for(int i=1; i<=pop_size; i++)
            for(int j=0; j<L; j++) {
                pro=rand()%1000/1000.0; //�ڡ�0,1��������������
                if(pro<pm) {
                    bitset<L>v(p[i]);
                    v.flip(j);
                    p[i]=v.to_ulong();
                }
            }
        t++;
        cout<<endl; //��Ⱥ��ֳһ��
    }
    cout<<"���ս����";
    for(int i(0); i<pop_size; i++) { //�㷨������������
        cout<<p[i]<<" ";
    }
    cout<<endl;
    return 0;
}
//�����ӽ�����
void Xover(int &a,int &b) {
    int pos; //��������ӽ��� ���ڼ������������໥����
    pos=rand()%5+1; //��n�������У����ȷ����pos������
    int j,k;
    j=pos;
    k=pos;
    bitset<L>e(a);
    bitset<L>f(b); //ǰpos�����������໥����
    bitset<L>g;
    bitset<L>h;
    for(int i=0; i<pos; i++) {
        if(e[i]==1)
            g.set(i);
    }
    for(int i=0; i<pos; i++) {
        if(f[i]==1)
            h.set(i);
    }
    for(j; j<L; j++) {
        if(f[j]==1)
            g.set(j);
    }
    for(k; k<L; k++) {
        if(e[k]==1)
            h.set(k);
    }
    a=g.to_ulong();
    b=h.to_ulong();
}
