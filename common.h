#include <memory>
struct free_deleter{
    template <typename T>
    void operator()(T *p) const {
        std::free(const_cast<std::remove_const_t<T>*>(p));
    }
};
template <typename T>
class unique_C_ptr: public std::unique_ptr<T[],free_deleter>{
  public:
    using std::unique_ptr<T[],free_deleter>::unique_ptr;
    // operator T*(){
    //     return c_ptr();

    // }
    void realloc(size_t new_size){
        auto ptr = this->release();
        ptr = (T*)::realloc(ptr, new_size);
        this->reset(ptr);
    }
};

static_assert(sizeof(char *)==
              sizeof(unique_C_ptr<char>),""); // ensure no overhead

template<typename T> 
auto make_malloc(size_t size) -> unique_C_ptr<T>{
    return unique_C_ptr<T>((T*)malloc(size * sizeof(T)));
}

template<typename T> 
auto make_calloc(size_t size) -> unique_C_ptr<T>{
    return unique_C_ptr<T>((T*)calloc(size, sizeof(T)));
}

