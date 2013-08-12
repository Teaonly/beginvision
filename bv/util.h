#ifndef _BV_UTIL_H_
#define _BV_UTIL_H_

namespace bv {

template<typename T> 
T min(T& v1, T& v2) {
    if ( v1 < v2)
        return v1;

    return v2;
}


}
#endif
