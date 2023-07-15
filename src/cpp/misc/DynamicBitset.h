//
// Created by knirschl on 06.07.23.
//

#ifndef BA_DYNAMICBITSET_H
#define BA_DYNAMICBITSET_H

#include <vector>
#include <stdint.h>

/**
 * Implementation of a dynamic bitset by ASTER.
 * @link https://github.com/chaoszhang/ASTER/blob/master/src/astral-pro.cpp#L55
 */
struct DynamicBitset {
private:
    int size{};
    std::vector<uint64_t> vec;

public:
    DynamicBitset() = default;

    explicit DynamicBitset(int sz) : size(sz), vec((sz + 63) / 64) {}

    /**
     * Copy constructor.
     *
     * @param copy_db
     */
    DynamicBitset(const DynamicBitset &copy_db) : size{copy_db.size}, vec(size) {
        std::copy(copy_db.vec.begin(), copy_db.vec.end(), this->vec.begin());
    }

    /**
     * Set the i-th bit.
     *
     * @param i position of the bit to set
     */
    void set(int i) {
        if (i >= size) {
            size = i + 1;
            if ((size + 63) / 64 > vec.size()) {
                vec.resize((size + 63) / 64);
            }
        }
        vec[i / 64] |= (1LL << (i % 64));
    }

    /**
     * ORs two bitsets without changing them. The smaller bitset gets virtually padded with zeros at
     * the higher bits.
     *
     * @param b the other bitset
     * @return a new DynamicBitset of the size of the bigger one
     */
    DynamicBitset operator|(const DynamicBitset &b) const {
        if (size < b.size) {
            return b | *this;
        }
        DynamicBitset res(size);
        for (int i = 0; i < b.vec.size(); i++) {
            res.vec[i] = vec[i] | b.vec[i];
        }
        for (int i = b.vec.size(); i < vec.size(); i++) {
            res.vec[i] = vec[i];
        }
        return res;
    }

    /**
     * ANDs two bitsets without changing them. The bigger bitset gets virtually truncated at the
     * higher bits.
     *
     * @param b the other bitset
     * @return a new DynamicBitset of the size of the smaller one
     */
    DynamicBitset operator&(const DynamicBitset &b) const {
        if (size < b.size) {
            return b & *this;
        }
        DynamicBitset res(b.size);
        for (int i = 0; i < b.vec.size(); i++) {
            res.vec[i] = vec[i] & b.vec[i];
        }
        return res;
    }

    /**
     * XORs two bitsets without changing them. All bits at greater positions than the size of the
     * smaller bitset are copied from the bigger bitset.
     *
     * @param b the other bitset
     * @return a new DynamicBitset of the size of the bigger one
     */
    DynamicBitset operator^(const DynamicBitset &b) const {
        if (size < b.size) {
            return b ^ *this;
        }
        DynamicBitset res(size);
        for (int i = 0; i < b.vec.size(); i++) {
            res.vec[i] = vec[i] ^ b.vec[i];
        }
        for (int i = b.vec.size(); i < vec.size(); i++) {
            res.vec[i] = vec[i];
        }
        return res;
    }

    /**
     * Subtracts two bitsets without changing them. Subtracting two bits from one another is defined
     * as ANDing with the NEGated subtrahend.
     * @code a - b := a & ~b
     *
     * @param b the other bitset
     * @return a new DynamicBitset of the size of the minuend (this)
     */
    DynamicBitset operator-(const DynamicBitset &b) const {
        DynamicBitset res(size);
        for (int i = 0; i < vec.size(); i++) {
            if (i < b.vec.size()) {
                res.vec[i] = vec[i] & ~b.vec[i];
            } else {
                res.vec[i] = vec[i];
            }
        }
        return res;
    }

    /**
     * NEGates a bitset without changing it.
     *
     * @return a new DynamicBitset of the size of this bitset
     */
    DynamicBitset operator~() const {
        DynamicBitset res(size);
        for (int i = 0; i < vec.size(); i++) {
            res.vec[i] = ~vec[i];
        }
        return res;
    }

    /**
     * Compares two bitsets for equality. Two bitsets are equal if they have the same bit at all
     * positions and the bigger bitset consists of only zeros at every position greater than the
     * size of the smaller bitset.
     *
     * @param b the other bitset
     * @return true if the two DynamicBitsets are equal under the above definition, false otherwise
     */
    bool operator==(const DynamicBitset &b) const {
        if (size < b.size) {
            return b == *this;
        }
        for (int i = 0; i < b.vec.size(); i++) {
            if (vec[i] != b.vec[i]) {
                return false;
            }
        }
        for (int i = b.vec.size(); i < vec.size(); i++) {
            if (vec[i] != 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * Compares two bitsets for inequality. Two bitsets are inequal if they have a different bit at
     * all positions and the bigger bitset consists of only ones at every position greater than the
     * size of the smaller bitset.
     *
     * @param b the other bitset
     * @return true if the two DynamicBitsets are not equal under the above definition, false
     *          otherwise
     */
    bool operator!=(const DynamicBitset &b) const {
        return !(*this == b);
    }

    /**
     * Checks if two bitsets are disjoint to one another. Two bitsets are disjoint if at every
     * position only one or none of the bitsets has a bit set.
     *
     * @param b the other bitset
     * @return true if the two DynamicBitsets are disjoint to one another under the above
     *          definition, false otherwise
     */
    bool isDisjointTo(const DynamicBitset &b) const {
        if (size < b.size) {
            return b.isDisjointTo(*this);
        }
        for (int i = 0; i < b.vec.size(); i++) {
            if ((vec[i] & b.vec[i]) != 0) {
                return false;
            }
        }
        return true;
    }
};

#endif //BA_DYNAMICBITSET_H
