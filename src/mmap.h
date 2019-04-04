#include <vector>
#include <iostream>
#include <iomanip>
#include <sys/mman.h>
#include <fcntl.h>

const size_t page_size = 2UL * 1024 * 1024;

template<typename T>
class mmap_allocator {
public:
    typedef T value_type;
    T *allocate(size_t n) {
        size_t length = n * sizeof(T);
        if(allocated >= length){
            return reinterpret_cast<T *>(mm_data);
        }
        if(length % page_size != 0){
            length = (length / page_size + 1) * page_size;
        }
        allocated = length;
        mm_data = mmap(NULL, length, (PROT_READ | PROT_WRITE),
                       MAP_PRIVATE | MAP_ANONYMOUS | 0x40000, -1, 0);
        if(mm_data == NULL){
            perror("mmap");
        }
        return reinterpret_cast<T *>(mm_data);
    }

    void deallocate(T *p, size_t n) {
        if(mm_data != p){
            munmap(p, n);
        }
    }

private:
    void * mm_data;
    size_t allocated;
};

template<typename T>
class mm_vector : public std::vector<T, mmap_allocator<T>> {};
