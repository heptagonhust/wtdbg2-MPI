#include<iostream>
#include<fstream>
#include<cstring>
#include<vector>
#include<algorithm>

int main(int argc,char **argv){
    if(argc == 1){
        printf("Please enter the input file\n");
        return 0;
    }
    std::ifstream infile;
    infile.open(argv[1]);
    if(!infile.is_open()){
        printf("No such file\n");
        return 0;
    }
    int basecnt = 0;
    std::string in;
    std::vector<int>basenum;
    while(getline(infile,in)){
        if(in[0]=='>'){
            int addr = 0;
            for(auto x = in.begin();;){
                if(*x =='='){
                    in.erase(x);
                    break;
                }
                in.erase(x);
            }
            //std::cout<<in<<std::endl;
            basecnt += std::stoi(in);
            basenum.push_back(std::stoi(in));
        }else continue;
    }
    std::sort(basenum.begin(),basenum.end(),std::greater<int>());
    int nownum = 0;
    int N50 = 0;
    for(auto x:basenum){
        nownum += x;
        if(!N50 && x>=basecnt/2){
            printf("N50 = %d\n",x);
            N50 = 1;
        }
        if(x>=basecnt*0.9){
            printf("N90 = %d\n",x);
            break;
        }
    }
}