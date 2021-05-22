#include <iostream>
#include <cstring>
#include <random>
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <cassert>

#define Times 100
#define DATASEGMENTSIZE 512

using namespace std;

typedef struct Record{
    int64_t *seg_start=NULL;
    int records=0;
}record;

int64_t nextData(int64_t, double);
int64_t PMAWithBranching(int64_t *, int64_t, int64_t, int64_t);
int64_t PMAWithDirectBranching(int64_t *, int64_t, int64_t, int64_t);
int64_t PMAWithDenseData(int64_t *, int64_t, record *, int64_t, int64_t, int64_t);
int64_t PMAWithDenseData2(int64_t *, int64_t, record *, int64_t, int64_t, int64_t);
int64_t PMAWithPositionOffset(int64_t *, int64_t *, int64_t, int64_t, int64_t);
void create_next_element_offset();
void create_non_zero_entries();
int64_t PMAWithNextElementOffset(int64_t *, unsigned short *, int64_t, int64_t, int64_t);
int64_t PMAWithNonZeroEntries(int64_t *, unsigned short *, int64_t, int64_t, int64_t);

unsigned short NextElementOffset[65536][16];
u_char *NonZeroEntries[65536];

int main(int argc, char **argv){
    int64_t *data;
    int64_t *index;
    unsigned short *NEOffset;
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
    memset(data, 0, sizeof(data));

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
    //bitIndex = (uint64_t *) malloc (sizeof(uint64_t) * ceil((double)count/sizeof(uint64_t)));
    for(int64_t i=0, ind = 0; i<arraySize; i++){
        if(data[i] > 0){
            index[ind++] = i;
        }
    }

    cout <<"Creating bitmap and next element offset for data." << endl;
    create_next_element_offset();
    int ushortSize = sizeof(unsigned short) * 8;
    int64_t NEOffset_size = ceil((double)arraySize / ushortSize);
    NEOffset = (unsigned short *) malloc(sizeof(unsigned short) * NEOffset_size);
    for(int64_t i=0; i<NEOffset_size; i++){
        int64_t start = i * ushortSize;
        int64_t limit = min((int64_t)start+ushortSize, arraySize);
        unsigned short k = 0, l = 1;
        //if(i == 2095059) cout <<"start: "<<start<<" limit: " <<limit << " arraySize: "<<arraySize<<" ushortsize: "<<ushortSize<<" ";
        for(int64_t j=start; j<limit; j++){
            if(data[j] != 0) k = k|l;
            l = l << 1;
        }
        NEOffset[i] = k;
    }

    cout <<"Creating next nonzero element list for data." << endl;
    create_non_zero_entries();

    cout<< "creating dense array for Bconz" << endl;
    data2 = (int64_t *)malloc(sizeof(int64_t)*arraySize);
    memset(data2, 0, sizeof(data2));
    int64_t segments = ceil((double)arraySize/DATASEGMENTSIZE);
    record *segRecord = new record[segments];
    for(int64_t segNo=0; segNo < segments; segNo++){
        int64_t start = segNo * DATASEGMENTSIZE;
        int64_t limit = min(start + DATASEGMENTSIZE, arraySize);
        int count = 0;
        segRecord[segNo].seg_start = &data2[start];
        for(int64_t i = start, j = start; i<limit; i++){
            if(data[i] == 0) continue;
            data2[j++]= data[i];
            count++;
        }
        segRecord[segNo].records = count;
    }

    chrono::time_point<std::chrono::high_resolution_clock> start, end;
    cout << "Running "<<Times<<" short range (length 1000) queries." << endl;
    int64_t pamTimerShort = 0, bitOffTimerShort = 0, posOffTimerShort = 0, nzOffTimerShort = 0, brnchTimerShort = 0, dBrnchTimerShort = 0;
    for(int i = 0; i < Times; i++){
        //create range
        d = distribution(generator);
        int64_t low = (int64_t)(d*size) % current;
        if (low + 1000 > current) low -= 1000;

        //Run PMA with branching
        start = chrono::high_resolution_clock::now();
        int64_t valBR = PMAWithBranching(data, arraySize, low, low+1000);
        end = chrono::high_resolution_clock::now();
        //cout <<"Vaue at end "<<end<<endl;
        brnchTimerShort += chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        //Run PMA with direct access
        start = chrono::high_resolution_clock::now();
        valBR = PMAWithDirectBranching(data, arraySize, low, low+1000);
        end = chrono::high_resolution_clock::now();
        //cout <<"Vaue at end "<<end<<endl;
        dBrnchTimerShort += chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        //Run PMA with Boncz Dense Data
        start = chrono::high_resolution_clock::now();
        int64_t valB = PMAWithDenseData2(data2, arraySize, segRecord, segments, low, low+1000);
        end = chrono::high_resolution_clock::now();
        //cout <<"Vaue at end "<<end<<endl;
        pamTimerShort += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        assert(valB == valBR);

        //Run PMA with bitwise offset
        start = chrono::high_resolution_clock::now();
        int64_t valBO = PMAWithNextElementOffset(data, NEOffset, arraySize, low, low+1000);
        end = chrono::high_resolution_clock::now();
        bitOffTimerShort += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //printf("value from NEO %ld, value from Dense Data: %ld, value from branching: %ld\n",valBO, valB, valBR);
        assert(valBO == valB);

        //Run PMA with position offset
        start = chrono::high_resolution_clock::now();
        int64_t valPO = PMAWithPositionOffset(data, index, count, low, low+1000);
        end = chrono::high_resolution_clock::now();
        posOffTimerShort += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //printf("value from NEO %ld, value from Dense Data: %ld, value from Position offset: %ld\n",valBO, valB, valPO);
        assert(valBO == valPO);

        //Run PMA with Non-zero Entries
        start = chrono::high_resolution_clock::now();
        int64_t valNZO = PMAWithNonZeroEntries(data, NEOffset, arraySize, low, low+1000);
        end = chrono::high_resolution_clock::now();
        nzOffTimerShort += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //printf("value from NEO %ld, value from Dense Data: %ld, value from Position offset: %ld, value from Nonzero offset %ld\n",valBO, valB, valPO, valNZO);
        //assert(valBO == valNZO);
    }

    cout << "Running "<<Times<<" long range (length 100000) queries." << endl;
    int64_t pamTimerLong = 0, posOffTimerLong = 0, bitOffTimerLong = 0, nzOffTimerLong = 0, brnchTimerLong = 0, dBrnchTimerLong = 0;
    for(int i = 0; i < Times; i++){
        //create range
        d = distribution(generator);
        int64_t low = (int64_t)(d*size) % current;
        if (low + 100000 > current) low -= 100000;

        //Run PAM with branching
        start = chrono::high_resolution_clock::now();
        int64_t valBR = PMAWithBranching(data, arraySize, low, low+100000);
        end = chrono::high_resolution_clock::now();
        brnchTimerLong += chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        //Run PAM with direct access
        start = chrono::high_resolution_clock::now();
        valBR = PMAWithDirectBranching(data, arraySize, low, low+100000);
        end = chrono::high_resolution_clock::now();
        dBrnchTimerLong += chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        //Run PAM with Boncz Dense Data
        start = chrono::high_resolution_clock::now();
        int64_t valB = PMAWithDenseData2(data2, arraySize, segRecord, segments, low, low+100000);
        end = chrono::high_resolution_clock::now();
        pamTimerLong += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        assert(valB == valBR);

        //Run PAM with Bitwise offset
        start = chrono::high_resolution_clock::now();
        int64_t valBO = PMAWithNextElementOffset(data, NEOffset, arraySize, low, low+100000);
        end = chrono::high_resolution_clock::now();
        bitOffTimerLong += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //printf("value from NEO %ld, value from Dense Data: %ld, value from branching: %ld\n",valBO, valB, valBR);
        assert(valBO == valB);

        //Run PAM with Position offset
        start = chrono::high_resolution_clock::now();
        int64_t valPO = PMAWithPositionOffset(data, index, count, low, low+100000);
        end = chrono::high_resolution_clock::now();
        posOffTimerLong += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        assert(valBO == valPO);

        //Run PMA with Non-zero Entries
        start = chrono::high_resolution_clock::now();
        int64_t valNZO = PMAWithNonZeroEntries(data, NEOffset, arraySize, low, low+100000);
        end = chrono::high_resolution_clock::now();
        nzOffTimerLong += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //assert(valBO == valNZO);
        //printf("value from NEO %ld, value from Dense Data: %ld, value from Position offset: %ld, value from Nonzero offset %ld\n",valBO, valB, valPO, valNZO);
    }

    cout << "Running "<<Times<<" point lookup queries." << endl;
    int64_t pamTimerPoint = 0, posOffTimerPoint = 0, bitOffTimerPoint = 0, nzOffTimerPoint = 0, brnchTimerPoint = 0;
    for(int i = 0; i < Times; i++){
        //create range
        d = distribution(generator);
        int64_t low = (int64_t)(d*size) % current;

        //Run PMA with branching
        start = chrono::high_resolution_clock::now();
        int64_t valBR = PMAWithBranching(data, arraySize, low, 0);
        end = chrono::high_resolution_clock::now();
        brnchTimerPoint += chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        //Run PMA with Boncz dense data
        start = chrono::high_resolution_clock::now();
        int64_t valB = PMAWithDenseData2(data2, arraySize, segRecord, segments, low, 0);
        end = chrono::high_resolution_clock::now();
        pamTimerPoint += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        assert(valB == valBR);

        //Run PMA with btiwise offset
        start = chrono::high_resolution_clock::now();
        int64_t valBO = PMAWithNextElementOffset(data, NEOffset, arraySize, low, 0);
        end = chrono::high_resolution_clock::now();
        bitOffTimerPoint += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //printf("value from NEO %ld, value from Dense Data: %ld, value from branching: %ld\n",valBO, valB, valBR);
        assert(valBO == valB);

        //Run PMA with position offset
        start = chrono::high_resolution_clock::now();
        int64_t valPO = PMAWithPositionOffset(data, index, count, low, 0);
        end = chrono::high_resolution_clock::now();
        posOffTimerPoint += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        assert(valBO == valPO);

        //Run PMA with Non-zero Entries
        start = chrono::high_resolution_clock::now();
        int64_t valNZO = PMAWithNonZeroEntries(data, NEOffset, arraySize, low, 0);
        end = chrono::high_resolution_clock::now();
        nzOffTimerPoint += chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //printf("value from NEO %ld, value from Dense Data: %ld, value from Position offset: %ld, value from Nonzero offset %ld\n",valBO, valB, valPO, valNZO);
        assert(valBO == valNZO);
    }
    
    cout << "Average time for short range: PMA with Dense Data (Boncz) = "<<pamTimerShort<<", PMA with Bitmap offset ="<<bitOffTimerShort<<", PMA with Position offset ="<<posOffTimerShort<<", PMA with non-zero offset ="<<nzOffTimerShort<<endl;
    cout << "Average time for long range: PMA with Dense Data (Boncz) = "<<pamTimerLong<<", PMA with Bitmap offset ="<<bitOffTimerLong<<", PMA with Position offset ="<<posOffTimerLong<<", PMA with non-zero offset ="<<nzOffTimerLong<<endl;
    cout << "Average time for point lookup: PMA with Dense Data (Boncz) = "<<pamTimerPoint<<", PMA with Bitmap offset ="<<bitOffTimerPoint<<", PMA with Position offset ="<<posOffTimerPoint<<", PMA with non-zero offset ="<<nzOffTimerPoint<<endl;

    ofstream myfile;
    myfile.open ("Compare_short_range.csv", ios::out | ios::app);
    myfile << "Element," << size<< ",Density,"<<density<<",Repeated,"<<Times<<",Branching,"<<brnchTimerShort<<",DirectAccess,"<<dBrnchTimerShort<<",Boncz,"<<pamTimerShort<<",BitOffset,"<<bitOffTimerShort<<",PosOffst,"<<posOffTimerShort<<",NonZeroOffst,"<<nzOffTimerShort<<endl;
    myfile.close();
    myfile.open ("Compare_long_range.csv", ios::out | ios::app);
    myfile << "Element," << size<< ",Density,"<<density<<",Repeated,"<<Times<<",Branching,"<<brnchTimerLong<<",DirectAccess,"<<dBrnchTimerLong<<",Boncz,"<<pamTimerLong<<",BitOffst,"<<bitOffTimerLong<<",PosOffst,"<<posOffTimerLong<<",NonZeroOffst,"<<nzOffTimerLong<<endl;
    myfile.close();
    myfile.open ("Compare_lookup.csv", ios::out | ios::app);
    myfile << "Element," << size<< ",Density,"<<density<<",Repeated,"<<Times<<",Branching,"<<brnchTimerPoint<<",Boncz,"<<pamTimerPoint<<",BitOffst,"<<bitOffTimerPoint<<",PosOffst,"<<posOffTimerPoint<<",NonZeroOffst,"<<nzOffTimerPoint<<endl;
    myfile.close();

    delete[] data;  
    delete[] data2;
    delete[] index;
    delete[] NEOffset;
    for(int i=0;i<65536; i++) delete[] NonZeroEntries[i];
    return 0;
}

int64_t nextData(int64_t current, double d){
    return (int64_t)(d*(1-d)*20) + current + 1;
}

int64_t PMAWithBranching(int64_t *data, int64_t end, int64_t low, int64_t high){
    //find first element using binary search
    int64_t start = 0; end--;
    int64_t mid;
    //cout << "starting Branching with start "<<start<<" and end. "<<end<<"Searching for "<<low<<endl; 
    while(start <= end){
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

int64_t PMAWithDirectBranching(int64_t *data, int64_t end, int64_t low, int64_t high){
    //find first element using binary search
    int64_t start = 0; end--;
    int64_t mid;
    //cout << "starting Branching with start "<<start<<" and end. "<<end<<"Searching for "<<low<<endl; 
    while(start <= end){
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
        //if(data[mid] != 0) {
            //cout << "Index "<<mid<<" value"<<data[mid] <<" Sum: "<< sum<< endl;
            sum += data[mid];
        //}
        mid++;
    }
    //cout << "finished with "<<sum<<endl;
    return sum;
}

int64_t PMAWithPositionOffset(int64_t *data, int64_t *index, int64_t end, int64_t low, int64_t high){
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

int64_t PMAWithNextElementOffset(int64_t *data, unsigned short *bitMap, int64_t end, int64_t low, int64_t high){
    int64_t start = 0; end--;
    int64_t mid;
    while(start <= end){
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
    while(data[mid]<low){
        mid++;
    }
    if (high == 0) return data[mid]==low? low : 0;

    int ushortSize = sizeof(unsigned short) * 8;
    int64_t segNo = mid / ushortSize;
    unsigned short temp = bitMap[segNo];
    int64_t base;
    int64_t pBase = segNo*ushortSize;
    int64_t sum = 0;
    int count = 0;

    //add elements from current segment
    bool flag = false;
    for(int offset = 0; offset<ushortSize; offset += NextElementOffset[temp][offset]){
        base = pBase+offset;
        if(flag && data[base]> high) flag = false;
        else if(!flag && data[base] >= low) flag = true;
        if(flag) {count++; sum += data[base];}
    }

    //add elements from next segments
    while(true){
        pBase += ushortSize;
        //base = pBase;
        segNo++;
        temp = bitMap[segNo];
        for(int offset = 0; offset<ushortSize; offset += NextElementOffset[temp][offset]){
            base = pBase + offset;
            if(data[base] > high) return sum;
            sum += data[base];
            count++;
        }
    }

    return sum;
}

int64_t PMAWithNonZeroEntries(int64_t *data, unsigned short *bitMap, int64_t end, int64_t low, int64_t high){
    int64_t start = 0; end--;
    int64_t mid;
    while(start <= end){
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
    while(data[mid]<low){
        mid++;
    }
    if (high == 0) return data[mid]==low? low : 0;

    int ushortSize = sizeof(unsigned short) * 8;
    int64_t segNo = mid / ushortSize;
    unsigned short temp = bitMap[segNo];
    int64_t base;
    int64_t pBase = segNo * ushortSize;
    int64_t sum = 0;
    int count = 0;

    //add elements from current segment
    bool flag = false;
    for(int offset = 0; offset < NonZeroEntries[temp][0] ; ){
        offset++;
        base = pBase + NonZeroEntries[temp][offset];
        if(flag && data[base]> high) flag = false;
        else if(!flag && data[base] >= low) flag = true;
        if(flag) {count++; sum += data[base];}
    }

    //add elements from next segments
    while(true){
        pBase += ushortSize;
        
        segNo++;
        //temp = bitMap[segNo];
        unsigned char * ar = NonZeroEntries[bitMap[segNo]];
        int64_t * br = &data[pBase];
        //unsigned char *cr = &ar[1];
        int offset = ar[0];
        switch (offset){
            case 1:
                //sum += data[pBase+ar[1]];
                //sum += data[pBase+*(ar++)];
                //sum += *(br + *cr);
                sum += *(br + ar[1]);
                break;
            case 2:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                break;
            case 3:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                break;
            case 4:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                break;
            case 5:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                break;
            case 6:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                break;
            case 7:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                sum += *(br + ar[7]);
                break;
            case 8:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                sum += *(br + ar[7]);
                sum += *(br + ar[8]);
                break;
            case 9:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                sum += *(br + ar[7]);
                sum += *(br + ar[8]);
                sum += *(br + ar[9]);
                break;
            case 10:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                sum += *(br + ar[7]);
                sum += *(br + ar[8]);
                sum += *(br + ar[9]);
                sum += *(br + ar[10]);
                break;
            case 11:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                sum += *(br + ar[7]);
                sum += *(br + ar[8]);
                sum += *(br + ar[9]);
                sum += *(br + ar[10]);
                sum += *(br + ar[11]);
                break;
            case 12:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                sum += *(br + ar[7]);
                sum += *(br + ar[8]);
                sum += *(br + ar[9]);
                sum += *(br + ar[10]);
                sum += *(br + ar[11]);
                sum += *(br + ar[12]);
                break;
            case 13:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                sum += *(br + ar[7]);
                sum += *(br + ar[8]);
                sum += *(br + ar[9]);
                sum += *(br + ar[10]);
                sum += *(br + ar[11]);
                sum += *(br + ar[12]);
                sum += *(br + ar[13]);
                break;
            case 14:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                sum += *(br + ar[7]);
                sum += *(br + ar[8]);
                sum += *(br + ar[9]);
                sum += *(br + ar[10]);
                sum += *(br + ar[11]);
                sum += *(br + ar[12]);
                sum += *(br + ar[13]);
                sum += *(br + ar[14]);
                break;
            case 15:
                sum += *(br + ar[1]);
                sum += *(br + ar[2]);
                sum += *(br + ar[3]);
                sum += *(br + ar[4]);
                sum += *(br + ar[5]);
                sum += *(br + ar[6]);
                sum += *(br + ar[7]);
                sum += *(br + ar[8]);
                sum += *(br + ar[9]);
                sum += *(br + ar[10]);
                sum += *(br + ar[11]);
                sum += *(br + ar[12]);
                sum += *(br + ar[13]);
                sum += *(br + ar[14]);
                sum += *(br + ar[15]);
                break;
            case 16:
                sum += *(br + 1);
                sum += *(br + 2);
                sum += *(br + 3);
                sum += *(br + 4);
                sum += *(br + 5);
                sum += *(br + 6);
                sum += *(br + 7);
                sum += *(br + 8);
                sum += *(br + 9);
                sum += *(br + 10);
                sum += *(br + 11);
                sum += *(br + 12);
                sum += *(br + 13);
                sum += *(br + 14);
                sum += *(br + 15);
                sum += *(br + 16);
        }
	    if(*(br + ar[offset]) > high){
            while(offset > 0 && *(br + ar[offset]) > high){
                sum -= *(br + ar[offset]);
                offset--;
            }
		    return sum;
	    }
    }
    return sum;
}

int64_t PMAWithDenseData(int64_t *data, int64_t end, record *segRecord, int64_t segments, int64_t low, int64_t high){
    int64_t start = 0; end--;
    int64_t mid;
    //cout << "starting Branching with start "<<start<<" and end. "<<end<<"Searching for "<<low<<endl; 
    while(start <= end){
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
    //cout << "mid value "<<data[mid]<<" mid position: "<<mid<<endl;
    int64_t sum = 0;
    while(data[mid]<low){
        mid++;
    }
    if (high == 0) return data[mid]==low? low : 0;

        //Find the locatio of the start of the segment
    int pos = mid;
    while(data[pos]!=0 && pos >= 0){
        pos--;
    }
    pos++;

    int segNo = 0;
    //Find the appropriate segment, using binary search
    start = 0, end = segments;
    while(start<=end){
        segNo = (int)(start + end)/2;
        int64_t value = *segRecord[segNo].seg_start;
        if(value == data[pos]) break;
        else if(value < data[pos]) start = segNo + 1;
        else end = segNo - 1;
    }
    
    int64_t limit = segRecord[segNo].records - (mid-pos);
    //printf("initial limit: %ld\n",limit);
    while(true){
        while(limit>0){
            if(data[mid]<=high) sum+= data[mid];
            else return sum;
            mid++; limit--;
        }
        segNo++;
        limit = segRecord[segNo].records;
        mid = pos + DATASEGMENTSIZE;
        pos = mid;
    }
    //cout << "finished with "<<sum<<endl;
    return sum;
}

int64_t PMAWithDenseData2(int64_t *data, int64_t end, record *segRecord, int64_t segments, int64_t low, int64_t high){
    int64_t segNo = 0, position = 0, segStart = 0, first = 0, last = segments - 1;
    int64_t sum = 0;
    while(first<=last){
        segNo = (first + last)/2;
        if(*segRecord[segNo].seg_start <= low && *(segRecord[segNo].seg_start+(segRecord[segNo].records-1)) >= low) break;
        else if(*segRecord[segNo].seg_start > low) last = segNo - 1;
        else first = segNo + 1;
    }
    //Do binary search later
    segStart = segNo * DATASEGMENTSIZE;
    first = 0, last = segRecord[segNo].records-1;
    
    for(int64_t i = 0; i < segRecord[segNo].records; i++){
        if(data[i+segStart] >= low) {position = i; break;}
    }
    
    if (high == 0) return data[segStart+position]==low? low : 0;

    int64_t limit = segRecord[segNo].records - position + 1;
    position += segStart;
    //printf("initial limit: %ld\n",limit);
    while(true){
        while(limit>0){
            if(data[position] > high) return sum;
            sum+= data[position];
            position++; limit--;
        }
        segNo++;
        limit = segRecord[segNo].records;
        position = segStart + DATASEGMENTSIZE;
        segStart = position;
    }
    //cout << "finished with "<<sum<<endl;
    return sum;
}

void create_next_element_offset(){
    for(int i = 0; i < 65536; i++){
        int j = 0, k = 1, l = 1, n = 1;
        while(n < 16){
            k = k << 1;
            if(i & k){
                for(int m = j; m<n; m++){
                    NextElementOffset[i][m]=l;
                    l--;
                }
                j = n;
            }
            n++;
            l++;
        }
        for(int m = j; m < n; m++){
            NextElementOffset[i][m] = l--;
        }
    }
}

void create_non_zero_entries(){
    for(int i = 0; i<65536; i++){
        int j = 0, k = 1, count = 0, idx = 0;
        for(; j<16; j++){
            if(i & k) count++;
            k = k<<1;
        }
        NonZeroEntries[i] = new u_char[count+1];
        NonZeroEntries[i][idx++] = (u_char)count;
        k = 1;
        for(j = 0; j<16; j++){
            if(i & k) NonZeroEntries[i][idx++] = j;
            k = k<<1;
        }
        /*
        for(j=2; j<=NonZeroEntries[i][0]; j++){
            NonZeroEntries[i][j] -= NonZeroEntries[i][j-1];
        }
        */
    }
}
