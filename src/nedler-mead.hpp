#ifndef NELDER_MEAD
#define NELDER_MEAD
//------------------------------------------------------------------------------------------------------------------------------------------------
#include <numeric>
#include <algorithm>
#include <vector>
#include <tuple>
//------------------------------------------------------------------------------------------------------------------------------------------------
namespace angstrom {
namespace math {
namespace internal{
//------------------------------------------------------------------------------------------------------------------------------------------------
template<int TCount, class TBaseType>
class Essence {
private:
    template<int ...S>
    struct seq{};

    template<int N, int ...S>
    struct gen: gen<N - 1 , N - 1, S...> {};
    template<int ...S>
    struct gen<0, S...> {
      typedef seq<S...> type;
    };

    template<int N, class ...TBase>
    struct genArgs: genArgs<N - 1, TBaseType, TBase...>{};
    template<class ...TBase>
    struct genArgs<0, TBase...>{
        typedef std::tuple<TBase...> tuple;
    };


    template<class BinaryOperation, int ...S>
    inline TBaseType _call(seq<S...>, BinaryOperation op) {
        return op(std::get<S>(data)...);
    }

    template<class ...TArgs>
    inline void _empty(TArgs... args){}
    template<int S, class TVal, class Binaryopration>
    inline TVal _apply_(TVal &start, Binaryopration op){
        op(start, std::get<S>(data));
        return start;
    }
    template<class TVal, class BinaryOperation, int ...S>
    inline TVal _apply(seq<S...> , TVal &start, BinaryOperation op) {
        _empty(_apply_<S>(start, op)...);
        return start;
    }

    template<int S, class TVal, class Binaryopration>
    inline TVal _bin_(TVal &start, Essence &essence, Binaryopration op){
        op(start, std::get<S>(data), std::get<S>(essence.data));
        return start;
    }
    template<class TVal, class BinaryOperation, int ...S>
    inline TVal _bin(seq<S...> , TVal &start, Essence &essence, BinaryOperation op) {
        _empty(_bin_<S>(start, essence, op)...);
        return start;
    }
public:
    typedef typename genArgs<TCount>::tuple tuple;
private:
    tuple data;
public:
    Essence() = delete;
    template<class ...Args>
    Essence(Args... args):
        data(args...)
    {}
    template<class BinaryOpration>
    inline TBaseType call(BinaryOpration op) {
        return _call(typename gen<TCount>::type(), op);
    }
    template<class TVal, class BinaryOperation>
    inline TVal apply(TVal start, BinaryOperation op) {
        return _apply(typename gen<TCount>::type(), start, op);
    }
    template<class TVal, class BinaryOperation>
    inline TVal bin(TVal start, Essence &essence, BinaryOperation op) {
        return _bin(typename gen<TCount>::type(), start, essence, op);
    }
    Essence<TCount, TBaseType> operator+(Essence<TCount, TBaseType> essence) {
        Essence<TCount, TBaseType> result = *this;
        result.bin(0, essence, [](int &i, TBaseType &first, TBaseType &second){
            first += second;
        });
        return result;
    }
    Essence<TCount, TBaseType> operator-(Essence<TCount, TBaseType> essence) {
        Essence<TCount, TBaseType> result = *this;
        result.bin(0, essence, [](int &i, TBaseType &first, TBaseType &second){
            first -= second;
        });
        return result;
    }
    Essence<TCount, TBaseType> operator*(TBaseType item) {
        Essence<TCount, TBaseType> result = *this;
        result.apply(0, [item](int accumulate, TBaseType &elemet){
            elemet *= item;
        });
        return result;
    }
    Essence<TCount, TBaseType> operator/(TBaseType item) {
        Essence<TCount, TBaseType> result = *this;
        result.apply(0, [item](int accumulate, TBaseType &elemet){
            elemet /= item;
        });
        return result;
    }
    tuple getTuple(){
        return data;
    }
};
//------------------------------------------------------------------------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------------------------------------------------------------------------
template<int TCount>
class NelderMead {
    typedef internal::Essence<TCount, double> NelderC;
public:
    NelderMead() = delete;
    template<class ...TArgs>
    NelderMead(TArgs... begCoeff): bValue(begCoeff...),
        reflection(1.0),
        compress(0.5),
        streach(2.0)
    {}
    NelderMead(NelderMead &) = delete;

    template<class BinaryOperation>
    typename NelderC::tuple minimizing(double delta, BinaryOperation op){
        //algoritm http://strelka.ftf2.tsu.ru/~sid/mopt/node8.html
        NelderC start = bValue;
        NelderC result = bValue;
        //Creatre simplex of N-dimensional space
        std::vector<std::pair<double, NelderC>> simplex = start.apply(std::vector<std::pair<double, NelderC>>(),
                                                                      [&start, op](std::vector<std::pair<double, NelderC>> &accoumulate, double &param){
            accoumulate.push_back(std::pair<double, NelderC>(start.call(op), start));
            param += 1.;
        });
        simplex.push_back(std::pair<double, NelderC>(start.call(op), start));
        while(1) {
            //STEP 1
            //Calculate y = f(Xi)
            std::sort(simplex.begin(), simplex.end(), [op](std::pair<double, NelderC> first, std::pair<double, NelderC> second){
                return first.first < second.first;
            });

            std::pair<double, NelderC> xs = simplex.back();
            std::pair<double, NelderC> xr = simplex.front();

            NelderC xmidle = std::accumulate(simplex.begin() + 1,
                                             simplex.end() - 1,
                                             simplex.front().second / double(TCount),
                                             [](NelderC accumul, std::pair<double, NelderC> element){
                    return accumul + (element.second / double(TCount));
            });
            //Stoping calculating
            NelderC _xcenter = (xmidle * double(TCount) + simplex.back().second) / double(TCount + 1);
            std::pair<double, NelderC> xcenter(_xcenter.call(op), _xcenter);
            if(std::accumulate(simplex.begin(), simplex.end(), true, [&xcenter, delta](bool flag, std::pair<double, NelderC> element){
                return flag && (delta > fabs(element.first - xcenter.first) && element.first >= xcenter.first);
            })) {
                result = _xcenter;
                break;
            }
            //STEP 2
            NelderC xreflec = xmidle + (xmidle - xs.second) * reflection;
            if(xr.first > xreflec.call(op)) {
                //STEP 3
                NelderC xone = xmidle + (xreflec - xmidle) * streach;
                simplex.back() = std::pair<double, NelderC>(xreflec.call(op), xreflec);
                if(xreflec.call(op) > xone.call(op))
                    simplex.back() = std::pair<double, NelderC>(xone.call(op), xone);;
                continue;
            }
            //STEP 4
            NelderC secondMax = (*(simplex.end() - 2)).second;
            if(secondMax.call(op) > xreflec.call(op)){
                simplex.back() = std::pair<double, NelderC>(xreflec.call(op), xreflec);
                continue;
            }
            //STEP 5
            NelderC xI = std::min(xreflec, xs.second, [op](NelderC first, NelderC second){
                return first.call(op) < second.call(op);
            });
            NelderC xII = xmidle + (xI - xmidle) * compress;
            if(xII.call(op) > xI.call(op))
                std::for_each(simplex.begin(), simplex.end(), [op, &xr](std::pair<double, NelderC> &pair){
                    NelderC element = pair.second;
                    element = element + (xr.second - element) * 0.5;
                    pair = std::pair<double, NelderC>(element.call(op), element);
                });
            else
                simplex.back() = std::pair<double, NelderC>(xII.call(op), xII);
        }
        return  result.getTuple();
    }
private:
    const double reflection, compress, streach;
    NelderC bValue;
};
//------------------------------------------------------------------------------------------------------------------------------------------------
}
}
#endif

