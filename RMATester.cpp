#include <iostream>
#include <cstring>
#include <random>
#include <chrono>
#include <unistd.h>
#include <fstream>

#define Times 100
#define DATASEGMENTSIZE 256

using namespace std;

int64_t nextData(int64_t, double);
int64_t PAMWithBranching(int64_t *, int64_t, int64_t, int64_t);
int64_t PAMWithDenseData(int64_t *, int64_t, int64_t, int64_t);
int64_t PAMWithPositionOffset(int64_t *, int64_t *, int64_t, int64_t, int64_t);

int main(int argc, char **argv){
    int64_t *data;
    int64_t *index;
    uint64_t *bitIndex;
    int64_t *data2;
    int64_t size = 0;
    double density = 0.0;

    //provide #of elements and density as parameter
    if(argc<2) {
        cout << "Provide # of elements and data density as parameter."<<endl;
        exit(1);
    }
    size = atoi(argv[1]);
    density = atof(argv[2]);
    int64_t arraySize = size/density;

    cout<< "Creating array of size "<< arraySize << " for "<< size << " elements with density "<< density<< endl;
    data = (int64_t *)malloc(sizeof(int64_t)*arraySize);
    memset(data, 0, sizeof(sizeof(int64_t)*arraySize));

    cout<< "Pushing data in the array"<< endl;
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);
    double d = 0.0;
    int64_t current = 1, count = 0;
    for(int64_t i = 0; i < arraySize; i++ ){
        d = distribution(generator);
        if(d > 1-density){
            data[i] = nextData(current, d);
            current = data[i];
            count++;
        }
    }

    cout <<"Creating index offset for data." << endl;
    index = (int64_t *) malloc (sizeof(int64_t)*(count+1));
    bitIndex = (uint64_t *) malloc (sizeof(uint64_t) * ceil((double)count/sizeof(uint64_t)));
    memset(bitIndex, 0, sizeof(uint64_t) * ceil((double)count/sizeof(uint64_t)));
    uint64_t bitind = 0;
    for(int64_t i=0, ind = 0; i<arraySize; i++){
        bitIndex[bitind] = bitIndex[bitind] << 1;
        if(data[i] > 0){
            index[ind++] = i;
            bitIndex[bitind] = bitIndex[bitind] | 1;
        }
        if((i+1) % 64 == 0) bitind++;
    }

    cout<< "creating dense array for Bconz";
    data2 = (int64_t *)malloc(sizeof(int64_t)*arraySize);
    memset(data2, 0, sizeof(sizeof(int64_t)*arraySize));
    int64_t segments = ceil((double)arraySize/DATASEGMENTSIZE);
    for(int64_t segNo=0; segNo < segments; segNo++){
        int64_t start = segNo * DATASEGMENTSIZE;
        int64_t limit = min(start + DATASEGMENTSIZE, arraySize);
        for(int64_t i = start, j = start; i<limit; i++){
            if(data[i] == 0) continue;
            data2[j++]= data[i];
        }
    }


    cout << "Run queries now within range 1 and "<< current << endl;

    chrono::time_point<std::chrono::high_resolution_clock> start, end;
    cout << "Running "<<Times<<" short range (length 1000) queries." << endl;
    int64_t pamTimerShort = 0, posOffTimerShort = 0;
    for(int i = 0; i < Times; i++){
        //create range
        d = distribution(generator);
        int64_t low = (int64_t)(d*size) % current;
        if (low + 1000 > current) low -= 1000;

        //Run PAM with branching
        start = chrono::high_resolution_clock::now();
        int64_t valB = PAMWithDenseData(data2, arraySize, low, low+100);
        end = chrono::high_resolution_clock::now();
        //cout <<"Vaue at end "<<end<<endl;
        pamTimerShort += chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        //Run PAM with positional offset
        start = chrono::high_resolution_clock::now();
        int64_t valPO = PAMWithPositionOffset(data, index, count, low, low+100);
        end = chrono::high_resolution_clock::now();
        posOffTimerShort += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    }

    cout << "Running "<<Times<<" long range (length 100000) queries." << endl;
    int64_t pamTimerLong = 0, posOffTimerLong = 0;
    for(int i = 0; i < Times; i++){
        //create range
        d = distribution(generator);
        int64_t low = (int64_t)(d*size) % current;
        if (low + 100000 > current) low -= 100000;

        //Run PAM with branching
        start = chrono::high_resolution_clock::now();

        int64_t valB = PAMWithDenseData(data2, arraySize, low, low+100000);
        end = chrono::high_resolution_clock::now();
        pamTimerLong += chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        //Run PAM with positional offset
        start = chrono::high_resolution_clock::now();
        int64_t valPO = PAMWithPositionOffset(data, index, count, low, low+100000);
        end = chrono::high_resolution_clock::now();
        posOffTimerLong += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    }

    cout << "Running "<<Times<<" point lookup queries." << endl;
    int64_t pamTimerPoint = 0, posOffTimerPoint = 0;
    for(int i = 0; i < Times; i++){
        //create range
        d = distribution(generator);
        int64_t low = (int64_t)(d*size) % current;

        //Run PAM with branching
        start = chrono::high_resolution_clock::now();
        int64_t valB = PAMWithDenseData(data2, arraySize, low, 0);
        end = chrono::high_resolution_clock::now();
        pamTimerPoint += chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        //Run PAM with positional offset
        start = chrono::high_resolution_clock::now();
        int64_t valPO = PAMWithPositionOffset(data, index, count, low, 0);
        end = chrono::high_resolution_clock::now();
        posOffTimerPoint += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    }
    cout << "Average time for short range: PAM with Branching = "<<pamTimerShort<<", PAM with index offset ="<<posOffTimerShort<<endl;
    cout << "Average time for long range: PAM with Branching = "<<pamTimerLong<<", PAM with index offset ="<<posOffTimerLong<<endl;
    cout << "Average time for point lookup: PAM with Branching = "<<pamTimerPoint<<", PAM with index offset ="<<posOffTimerPoint<<endl;
    ofstream myfile;
    myfile.open ("Boncz_short_range.csv", ios::out | ios::app);
    myfile << "Element," << size<< ",Density,"<<density<<",Repeated,"<<Times<<",P_Boncz,"<<pamTimerShort<<",P_Offst,"<<posOffTimerShort<<endl;
    myfile.close();
    myfile.open ("Boncz_long_range.csv", ios::out | ios::app);
    myfile << "Element," << size<< ",Density,"<<density<<",Repeated,"<<Times<<",P_Boncz,"<<pamTimerLong<<",P_Offst,"<<posOffTimerLong<<endl;
    myfile.close();
    myfile.open ("Boncz_lookup.csv", ios::out | ios::app);
    myfile << "Element," << size<< ",Density,"<<density<<",Repeated,"<<Times<<",P_Boncz,"<<pamTimerPoint<<",P_Offst,"<<posOffTimerPoint<<endl;
    myfile.close();
    return 0;
}

int64_t nextData(int64_t current, double d){
    return (int64_t)(d*(1-d)*20) + current + 1;
}

int64_t PAMWithBranching(int64_t *data, int64_t end, int64_t low, int64_t high){
    //find first element using binary search
    int64_t start = 0; end--;
    int64_t mid;
    //cout << "starting Branching with start "<<start<<" and end. "<<end<<"Searching for "<<low<<endl; 
    while(start < end){
        mid = (start + end) / 2;
        if(data[mid] == 0){
            int64_t changedMid = mid, offset = -1;
            while(changedMid > start){
                changedMid += offset;
                if(data[changedMid]!=0) break;
            }
            if(changedMid == start){
                if(data[start] >= low) {mid = start; break;}
                changedMid = mid;
                offset = 1;
                while(changedMid < end){
                    changedMid += offset;
                    if(data[changedMid]!=0) break;
                }
                if(changedMid == end){
                    mid = end; break;
                }else mid = changedMid;
            }else mid = changedMid;

        }
        if(data[mid] == low) break;
        else if(data[mid] < low) start = mid + 1;
        else end = mid - 1;
    }
    //cout << "mid value "<<data[mid]<<endl;
    int64_t sum = 0;
    while(data[mid]<low){
        mid++;
    }
    if (high == 0) return data[mid]==low? low : 0;
    while(data[mid] <= high){
        if(data[mid] != 0) {
            //cout << "Index "<<mid<<" value"<<data[mid] <<" Sum: "<< sum<< endl;
            sum += data[mid];
        }
        mid++;
    }
    //cout << "finished with "<<sum<<endl;
    return sum;
}

int64_t PAMWithPositionOffset(int64_t *data, int64_t *index, int64_t end, int64_t low, int64_t high){
    int64_t start = 0; end--;
    //cout << "starting Offset with start "<<start<<" and end "<<end<<" search for: "<<low<<endl; 
    int64_t mid, value;
    while(start <= end){
        mid = (start + end) / 2;
        value = data[index[mid]];
        if(value == low) break;
        else if(value < low) start = mid + 1;
        else end = mid - 1;
    }
    //cout << "index of mid: "<< index[mid] <<" data at mid: " <<data[index[mid]]<< " search value: "<<low<<endl;
    if(data[index[mid]] < low) mid++;
    if(high == 0) return data[index[mid]] == low? low : 0;
    int64_t sum = 0;
    //cout << "Data mid: " <<data[index[mid]] << endl;
    while(data[index[mid]] <= high){
        //cout << "Index "<<index[mid]<<" value"<<data[index[mid]] <<" Sum: "<< sum<< endl;
        sum += data[index[mid]];
        mid++;
    }
    //cout << "finished with "<<sum<<endl;
    return sum;
}

int64_t PAMWithDenseData(int64_t *data, int64_t end, int64_t low, int64_t high){
    int64_t start = 0; end--;
    int64_t mid;
    //cout << "starting Branching with start "<<start<<" and end. "<<end<<"Searching for "<<low<<endl; 
    while(start < end){
        mid = (start + end) / 2;
        if(data[mid] == 0){
            int64_t changedMid = mid, offset = -1;
            while(changedMid > start){
                changedMid += offset;
                if(data[changedMid]!=0) break;
            }
            if(changedMid == start){
                if(data[start] >= low) {mid = start; break;}
                changedMid = mid;
                offset = 1;
                while(changedMid < end){
                    changedMid += offset;
                    if(data[changedMid]!=0) break;
                }
                if(changedMid == end){
                    mid = end; break;
                }else mid = changedMid;
            }else mid = changedMid;

        }
        if(data[mid] == low) break;
        else if(data[mid] < low) start = mid + 1;
        else end = mid - 1;
    }
    //cout << "mid value "<<data[mid]<<endl;
    int64_t sum = 0;
    while(data[mid]<low){
        mid++;
    }
    if (high == 0) return data[mid]==low? low : 0;
    while(data[mid] <= high){
        if(data[mid] != 0) {
            //cout << "Index "<<mid<<" value"<<data[mid] <<" Sum: "<< sum<< endl;
            sum += data[mid];
        }
        mid++;
    }
    //cout << "finished with "<<sum<<endl;
    return sum;
}