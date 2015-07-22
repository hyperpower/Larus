/************************
 //  \file   BigInt.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   22 avr. 2015 
 ***********************/

#include "BigInt.h"

namespace Larus
{

// Memory management definitions have moved to the bottom of NumberlikeArray.hh.

// The templates used by these constructors and converters are at the bottom of
// BigUnsigned.hh.

BigUnsigned::BigUnsigned(unsigned long x)
{
	initFromPrimitive(x);
}
BigUnsigned::BigUnsigned(unsigned int x)
{
	initFromPrimitive(x);
}
BigUnsigned::BigUnsigned(unsigned short x)
{
	initFromPrimitive(x);
}
BigUnsigned::BigUnsigned(long x)
{
	initFromSignedPrimitive(x);
}
BigUnsigned::BigUnsigned(int x)
{
	initFromSignedPrimitive(x);
}
BigUnsigned::BigUnsigned(short x)
{
	initFromSignedPrimitive(x);
}

unsigned long BigUnsigned::toUnsignedLong() const
{
	return convertToPrimitive<unsigned long>();
}
unsigned int BigUnsigned::toUnsignedInt() const
{
	return convertToPrimitive<unsigned int>();
}
unsigned short BigUnsigned::toUnsignedShort() const
{
	return convertToPrimitive<unsigned short>();
}
long BigUnsigned::toLong() const
{
	return convertToSignedPrimitive<long>();
}
int BigUnsigned::toInt() const
{
	return convertToSignedPrimitive<int>();
}
short BigUnsigned::toShort() const
{
	return convertToSignedPrimitive<short>();
}

// BIT/BLOCK ACCESSORS

void BigUnsigned::setBlock(Index i, Blk newBlock)
{
	if (newBlock == 0) {
		if (i < len) {
			blk[i] = 0;
			zapLeadingZeros();
		}
		// If i >= len, no effect.
	} else {
		if (i >= len) {
			// The nonzero block extends the number.
			allocateAndCopy(i + 1);
			// Zero any added blocks that we aren't setting.
			for (Index j = len; j < i; j++)
				blk[j] = 0;
			len = i + 1;
		}
		blk[i] = newBlock;
	}
}

/* Evidently the compiler wants BigUnsigned:: on the return type because, at
 * that point, it hasn't yet parsed the BigUnsigned:: on the name to get the
 * proper scope. */
BigUnsigned::Index BigUnsigned::bitLength() const
{
	if (isZero())
		return 0;
	else {
		Blk leftmostBlock = getBlock(len - 1);
		Index leftmostBlockLen = 0;
		while (leftmostBlock != 0) {
			leftmostBlock >>= 1;
			leftmostBlockLen++;
		}
		return leftmostBlockLen + (len - 1) * N;
	}
}

void BigUnsigned::setBit(Index bi, bool newBit)
{
	Index blockI = bi / N;
	Blk block = getBlock(blockI), mask = Blk(1) << (bi % N);
	block = newBit ? (block | mask) : (block & ~mask);
	setBlock(blockI, block);
}

// COMPARISON
BigUnsigned::CmpRes BigUnsigned::compareTo(const BigUnsigned &x) const
{
	// A bigger length implies a bigger number.
	if (len < x.len)
		return less;
	else if (len > x.len)
		return greater;
	else {
		// Compare blocks one by one from left to right.
		Index i = len;
		while (i > 0) {
			i--;
			if (blk[i] == x.blk[i])
				continue;
			else if (blk[i] > x.blk[i])
				return greater;
			else
				return less;
		}
		// If no blocks differed, the numbers are equal.
		return equal;
	}
}

// COPY-LESS OPERATIONS

/*
 * On most calls to copy-less operations, it's safe to read the inputs little by
 * little and write the outputs little by little.  However, if one of the
 * inputs is coming from the same variable into which the output is to be
 * stored (an "aliased" call), we risk overwriting the input before we read it.
 * In this case, we first compute the result into a temporary BigUnsigned
 * variable and then copy it into the requested output variable *this.
 * Each put-here operation uses the DTRT_ALIASED macro (Do The Right Thing on
 * aliased calls) to generate code for this check.
 *
 * I adopted this approach on 2007.02.13 (see Assignment Operators in
 * BigUnsigned.hh).  Before then, put-here operations rejected aliased calls
 * with an exception.  I think doing the right thing is better.
 *
 * Some of the put-here operations can probably handle aliased calls safely
 * without the extra copy because (for example) they process blocks strictly
 * right-to-left.  At some point I might determine which ones don't need the
 * copy, but my reasoning would need to be verified very carefully.  For now
 * I'll leave in the copy.
 */
#define DTRT_ALIASED(cond, op) \
	if (cond) { \
		BigUnsigned tmpThis; \
		tmpThis.op; \
		*this = tmpThis; \
		return; \
	}

void BigUnsigned::add(const BigUnsigned &a, const BigUnsigned &b)
{
	DTRT_ALIASED(this == &a || this == &b, add(a, b));
	// If one argument is zero, copy the other.
	if (a.len == 0) {
		operator =(b);
		return;
	} else if (b.len == 0) {
		operator =(a);
		return;
	}
	// Some variables...
	// Carries in and out of an addition stage
	bool carryIn, carryOut;
	Blk temp;
	Index i;
	// a2 points to the longer input, b2 points to the shorter
	const BigUnsigned *a2, *b2;
	if (a.len >= b.len) {
		a2 = &a;
		b2 = &b;
	} else {
		a2 = &b;
		b2 = &a;
	}
	// Set prelimiary length and make room in this BigUnsigned
	len = a2->len + 1;
	allocate(len);
	// For each block index that is present in both inputs...
	for (i = 0, carryIn = false; i < b2->len; i++) {
		// Add input blocks
		temp = a2->blk[i] + b2->blk[i];
		// If a rollover occurred, the result is less than either input.
		// This test is used many times in the BigUnsigned code.
		carryOut = (temp < a2->blk[i]);
		// If a carry was input, handle it
		if (carryIn) {
			temp++;
			carryOut |= (temp == 0);
		}
		blk[i] = temp; // Save the addition result
		carryIn = carryOut; // Pass the carry along
	}
	// If there is a carry left over, increase blocks until
	// one does not roll over.
	for (; i < a2->len && carryIn; i++) {
		temp = a2->blk[i] + 1;
		carryIn = (temp == 0);
		blk[i] = temp;
	}
	// If the carry was resolved but the larger number
	// still has blocks, copy them over.
	for (; i < a2->len; i++)
		blk[i] = a2->blk[i];
	// Set the extra block if there's still a carry, decrease length otherwise
	if (carryIn)
		blk[i] = 1;
	else
		len--;
}

void BigUnsigned::subtract(const BigUnsigned &a, const BigUnsigned &b)
{
	DTRT_ALIASED(this == &a || this == &b, subtract(a, b));
	if (b.len == 0) {
		// If b is zero, copy a.
		operator =(a);
		return;
	} else if (a.len < b.len)
		// If a is shorter than b, the result is negative.
		throw "BigUnsigned::subtract: "
				"Negative result in unsigned calculation";
	// Some variables...
	bool borrowIn, borrowOut;
	Blk temp;
	Index i;
	// Set preliminary length and make room
	len = a.len;
	allocate(len);
	// For each block index that is present in both inputs...
	for (i = 0, borrowIn = false; i < b.len; i++) {
		temp = a.blk[i] - b.blk[i];
		// If a reverse rollover occurred,
		// the result is greater than the block from a.
		borrowOut = (temp > a.blk[i]);
		// Handle an incoming borrow
		if (borrowIn) {
			borrowOut |= (temp == 0);
			temp--;
		}
		blk[i] = temp; // Save the subtraction result
		borrowIn = borrowOut; // Pass the borrow along
	}
	// If there is a borrow left over, decrease blocks until
	// one does not reverse rollover.
	for (; i < a.len && borrowIn; i++) {
		borrowIn = (a.blk[i] == 0);
		blk[i] = a.blk[i] - 1;
	}
	/* If there's still a borrow, the result is negative.
	 * Throw an exception, but zero out this object so as to leave it in a
	 * predictable state. */
	if (borrowIn) {
		len = 0;
		throw "BigUnsigned::subtract: Negative result in unsigned calculation";
	} else
		// Copy over the rest of the blocks
		for (; i < a.len; i++)
			blk[i] = a.blk[i];
	// Zap leading zeros
	zapLeadingZeros();
}

/*
 * About the multiplication and division algorithms:
 *
 * I searched unsucessfully for fast C++ built-in operations like the `b_0'
 * and `c_0' Knuth describes in Section 4.3.1 of ``The Art of Computer
 * Programming'' (replace `place' by `Blk'):
 *
 *    ``b_0[:] multiplication of a one-place integer by another one-place
 *      integer, giving a two-place answer;
 *
 *    ``c_0[:] division of a two-place integer by a one-place integer,
 *      provided that the quotient is a one-place integer, and yielding
 *      also a one-place remainder.''
 *
 * I also missed his note that ``[b]y adjusting the word size, if
 * necessary, nearly all computers will have these three operations
 * available'', so I gave up on trying to use algorithms similar to his.
 * A future version of the library might include such algorithms; I
 * would welcome contributions from others for this.
 *
 * I eventually decided to use bit-shifting algorithms.  To multiply `a'
 * and `b', we zero out the result.  Then, for each `1' bit in `a', we
 * shift `b' left the appropriate amount and add it to the result.
 * Similarly, to divide `a' by `b', we shift `b' left varying amounts,
 * repeatedly trying to subtract it from `a'.  When we succeed, we note
 * the fact by setting a bit in the quotient.  While these algorithms
 * have the same O(n^2) time complexity as Knuth's, the ``constant factor''
 * is likely to be larger.
 *
 * Because I used these algorithms, which require single-block addition
 * and subtraction rather than single-block multiplication and division,
 * the innermost loops of all four routines are very similar.  Study one
 * of them and all will become clear.
 */

/*
 * This is a little inline function used by both the multiplication
 * routine and the division routine.
 *
 * `getShiftedBlock' returns the `x'th block of `num << y'.
 * `y' may be anything from 0 to N - 1, and `x' may be anything from
 * 0 to `num.len'.
 *
 * Two things contribute to this block:
 *
 * (1) The `N - y' low bits of `num.blk[x]', shifted `y' bits left.
 *
 * (2) The `y' high bits of `num.blk[x-1]', shifted `N - y' bits right.
 *
 * But we must be careful if `x == 0' or `x == num.len', in
 * which case we should use 0 instead of (2) or (1), respectively.
 *
 * If `y == 0', then (2) contributes 0, as it should.  However,
 * in some computer environments, for a reason I cannot understand,
 * `a >> b' means `a >> (b % N)'.  This means `num.blk[x-1] >> (N - y)'
 * will return `num.blk[x-1]' instead of the desired 0 when `y == 0';
 * the test `y == 0' handles this case specially.
 */
inline BigUnsigned::Blk getShiftedBlock(const BigUnsigned &num,
		BigUnsigned::Index x, unsigned int y)
{
	BigUnsigned::Blk part1 =
			(x == 0 || y == 0) ? 0 : (num.blk[x - 1] >> (BigUnsigned::N - y));
	BigUnsigned::Blk part2 = (x == num.len) ? 0 : (num.blk[x] << y);
	return part1 | part2;
}

void BigUnsigned::multiply(const BigUnsigned &a, const BigUnsigned &b)
{
	DTRT_ALIASED(this == &a || this == &b, multiply(a, b));
	// If either a or b is zero, set to zero.
	if (a.len == 0 || b.len == 0) {
		len = 0;
		return;
	}
	/*
	 * Overall method:
	 *
	 * Set this = 0.
	 * For each 1-bit of `a' (say the `i2'th bit of block `i'):
	 *    Add `b << (i blocks and i2 bits)' to *this.
	 */
	// Variables for the calculation
	Index i, j, k;
	unsigned int i2;
	Blk temp;
	bool carryIn, carryOut;
	// Set preliminary length and make room
	len = a.len + b.len;
	allocate(len);
	// Zero out this object
	for (i = 0; i < len; i++)
		blk[i] = 0;
	// For each block of the first number...
	for (i = 0; i < a.len; i++) {
		// For each 1-bit of that block...
		for (i2 = 0; i2 < N; i2++) {
			if ((a.blk[i] & (Blk(1) << i2)) == 0)
				continue;
			/*
			 * Add b to this, shifted left i blocks and i2 bits.
			 * j is the index in b, and k = i + j is the index in this.
			 *
			 * `getShiftedBlock', a short inline function defined above,
			 * is now used for the bit handling.  It replaces the more
			 * complex `bHigh' code, in which each run of the loop dealt
			 * immediately with the low bits and saved the high bits to
			 * be picked up next time.  The last run of the loop used to
			 * leave leftover high bits, which were handled separately.
			 * Instead, this loop runs an additional time with j == b.len.
			 * These changes were made on 2005.01.11.
			 */
			for (j = 0, k = i, carryIn = false; j <= b.len; j++, k++) {
				/*
				 * The body of this loop is very similar to the body of the first loop
				 * in `add', except that this loop does a `+=' instead of a `+'.
				 */
				temp = blk[k] + getShiftedBlock(b, j, i2);
				carryOut = (temp < blk[k]);
				if (carryIn) {
					temp++;
					carryOut |= (temp == 0);
				}
				blk[k] = temp;
				carryIn = carryOut;
			}
			// No more extra iteration to deal with `bHigh'.
			// Roll-over a carry as necessary.
			for (; carryIn; k++) {
				blk[k]++;
				carryIn = (blk[k] == 0);
			}
		}
	}
	// Zap possible leading zero
	if (blk[len - 1] == 0)
		len--;
}

/*
 * DIVISION WITH REMAINDER
 * This monstrous function mods *this by the given divisor b while storing the
 * quotient in the given object q; at the end, *this contains the remainder.
 * The seemingly bizarre pattern of inputs and outputs was chosen so that the
 * function copies as little as possible (since it is implemented by repeated
 * subtraction of multiples of b from *this).
 *
 * "modWithQuotient" might be a better name for this function, but I would
 * rather not change the name now.
 */
void BigUnsigned::divideWithRemainder(const BigUnsigned &b, BigUnsigned &q)
{
	/* Defending against aliased calls is more complex than usual because we
	 * are writing to both *this and q.
	 *
	 * It would be silly to try to write quotient and remainder to the
	 * same variable.  Rule that out right away. */
	if (this == &q)
		throw "BigUnsigned::divideWithRemainder: Cannot write quotient and remainder into the same variable";
	/* Now *this and q are separate, so the only concern is that b might be
	 * aliased to one of them.  If so, use a temporary copy of b. */
	if (this == &b || &q == &b) {
		BigUnsigned tmpB(b);
		divideWithRemainder(tmpB, q);
		return;
	}

	/*
	 * Knuth's definition of mod (which this function uses) is somewhat
	 * different from the C++ definition of % in case of division by 0.
	 *
	 * We let a / 0 == 0 (it doesn't matter much) and a % 0 == a, no
	 * exceptions thrown.  This allows us to preserve both Knuth's demand
	 * that a mod 0 == a and the useful property that
	 * (a / b) * b + (a % b) == a.
	 */
	if (b.len == 0) {
		q.len = 0;
		return;
	}

	/*
	 * If *this.len < b.len, then *this < b, and we can be sure that b doesn't go into
	 * *this at all.  The quotient is 0 and *this is already the remainder (so leave it alone).
	 */
	if (len < b.len) {
		q.len = 0;
		return;
	}

	// At this point we know (*this).len >= b.len > 0.  (Whew!)

	/*
	 * Overall method:
	 *
	 * For each appropriate i and i2, decreasing:
	 *    Subtract (b << (i blocks and i2 bits)) from *this, storing the
	 *      result in subtractBuf.
	 *    If the subtraction succeeds with a nonnegative result:
	 *        Turn on bit i2 of block i of the quotient q.
	 *        Copy subtractBuf back into *this.
	 *    Otherwise bit i2 of block i remains off, and *this is unchanged.
	 *
	 * Eventually q will contain the entire quotient, and *this will
	 * be left with the remainder.
	 *
	 * subtractBuf[x] corresponds to blk[x], not blk[x+i], since 2005.01.11.
	 * But on a single iteration, we don't touch the i lowest blocks of blk
	 * (and don't use those of subtractBuf) because these blocks are
	 * unaffected by the subtraction: we are subtracting
	 * (b << (i blocks and i2 bits)), which ends in at least `i' zero
	 * blocks. */
	// Variables for the calculation
	Index i, j, k;
	unsigned int i2;
	Blk temp;
	bool borrowIn, borrowOut;

	/*
	 * Make sure we have an extra zero block just past the value.
	 *
	 * When we attempt a subtraction, we might shift `b' so
	 * its first block begins a few bits left of the dividend,
	 * and then we'll try to compare these extra bits with
	 * a nonexistent block to the left of the dividend.  The
	 * extra zero block ensures sensible behavior; we need
	 * an extra block in `subtractBuf' for exactly the same reason.
	 */
	Index origLen = len; // Save real length.
	/* To avoid an out-of-bounds access in case of reallocation, allocate
	 * first and then increment the logical length. */
	allocateAndCopy(len + 1);
	len++;
	blk[origLen] = 0; // Zero the added block.

	// subtractBuf holds part of the result of a subtraction; see above.
	Blk *subtractBuf = new Blk[len];

	// Set preliminary length for quotient and make room
	q.len = origLen - b.len + 1;
	q.allocate(q.len);
	// Zero out the quotient
	for (i = 0; i < q.len; i++)
		q.blk[i] = 0;

	// For each possible left-shift of b in blocks...
	i = q.len;
	while (i > 0) {
		i--;
		// For each possible left-shift of b in bits...
		// (Remember, N is the number of bits in a Blk.)
		q.blk[i] = 0;
		i2 = N;
		while (i2 > 0) {
			i2--;
			/*
			 * Subtract b, shifted left i blocks and i2 bits, from *this,
			 * and store the answer in subtractBuf.  In the for loop, `k == i + j'.
			 *
			 * Compare this to the middle section of `multiply'.  They
			 * are in many ways analogous.  See especially the discussion
			 * of `getShiftedBlock'.
			 */
			for (j = 0, k = i, borrowIn = false; j <= b.len; j++, k++) {
				temp = blk[k] - getShiftedBlock(b, j, i2);
				borrowOut = (temp > blk[k]);
				if (borrowIn) {
					borrowOut |= (temp == 0);
					temp--;
				}
				// Since 2005.01.11, indices of `subtractBuf' directly match those of `blk', so use `k'.
				subtractBuf[k] = temp;
				borrowIn = borrowOut;
			}
			// No more extra iteration to deal with `bHigh'.
			// Roll-over a borrow as necessary.
			for (; k < origLen && borrowIn; k++) {
				borrowIn = (blk[k] == 0);
				subtractBuf[k] = blk[k] - 1;
			}
			/*
			 * If the subtraction was performed successfully (!borrowIn),
			 * set bit i2 in block i of the quotient.
			 *
			 * Then, copy the portion of subtractBuf filled by the subtraction
			 * back to *this.  This portion starts with block i and ends--
			 * where?  Not necessarily at block `i + b.len'!  Well, we
			 * increased k every time we saved a block into subtractBuf, so
			 * the region of subtractBuf we copy is just [i, k).
			 */
			if (!borrowIn) {
				q.blk[i] |= (Blk(1) << i2);
				while (k > i) {
					k--;
					blk[k] = subtractBuf[k];
				}
			}
		}
	}
	// Zap possible leading zero in quotient
	if (q.blk[q.len - 1] == 0)
		q.len--;
	// Zap any/all leading zeros in remainder
	zapLeadingZeros();
	// Deallocate subtractBuf.
	// (Thanks to Brad Spencer for noticing my accidental omission of this!)
	delete[] subtractBuf;
}

/* BITWISE OPERATORS
 * These are straightforward blockwise operations except that they differ in
 * the output length and the necessity of zapLeadingZeros. */

void BigUnsigned::bitAnd(const BigUnsigned &a, const BigUnsigned &b)
{
	DTRT_ALIASED(this == &a || this == &b, bitAnd(a, b));
	// The bitwise & can't be longer than either operand.
	len = (a.len >= b.len) ? b.len : a.len;
	allocate(len);
	Index i;
	for (i = 0; i < len; i++)
		blk[i] = a.blk[i] & b.blk[i];
	zapLeadingZeros();
}

void BigUnsigned::bitOr(const BigUnsigned &a, const BigUnsigned &b)
{
	DTRT_ALIASED(this == &a || this == &b, bitOr(a, b));
	Index i;
	const BigUnsigned *a2, *b2;
	if (a.len >= b.len) {
		a2 = &a;
		b2 = &b;
	} else {
		a2 = &b;
		b2 = &a;
	}
	allocate(a2->len);
	for (i = 0; i < b2->len; i++)
		blk[i] = a2->blk[i] | b2->blk[i];
	for (; i < a2->len; i++)
		blk[i] = a2->blk[i];
	len = a2->len;
	// Doesn't need zapLeadingZeros.
}

void BigUnsigned::bitXor(const BigUnsigned &a, const BigUnsigned &b)
{
	DTRT_ALIASED(this == &a || this == &b, bitXor(a, b));
	Index i;
	const BigUnsigned *a2, *b2;
	if (a.len >= b.len) {
		a2 = &a;
		b2 = &b;
	} else {
		a2 = &b;
		b2 = &a;
	}
	allocate(a2->len);
	for (i = 0; i < b2->len; i++)
		blk[i] = a2->blk[i] ^ b2->blk[i];
	for (; i < a2->len; i++)
		blk[i] = a2->blk[i];
	len = a2->len;
	zapLeadingZeros();
}

void BigUnsigned::bitShiftLeft(const BigUnsigned &a, int b)
{
	DTRT_ALIASED(this == &a, bitShiftLeft(a, b));
	if (b < 0) {
		if (b << 1 == 0)
			throw "BigUnsigned::bitShiftLeft: "
					"Pathological shift amount not implemented";
		else {
			bitShiftRight(a, -b);
			return;
		}
	}
	Index shiftBlocks = b / N;
	unsigned int shiftBits = b % N;
	// + 1: room for high bits nudged left into another block
	len = a.len + shiftBlocks + 1;
	allocate(len);
	Index i, j;
	for (i = 0; i < shiftBlocks; i++)
		blk[i] = 0;
	for (j = 0, i = shiftBlocks; j <= a.len; j++, i++)
		blk[i] = getShiftedBlock(a, j, shiftBits);
	// Zap possible leading zero
	if (blk[len - 1] == 0)
		len--;
}

void BigUnsigned::bitShiftRight(const BigUnsigned &a, int b)
{
	DTRT_ALIASED(this == &a, bitShiftRight(a, b));
	if (b < 0) {
		if (b << 1 == 0)
			throw "BigUnsigned::bitShiftRight: "
					"Pathological shift amount not implemented";
		else {
			bitShiftLeft(a, -b);
			return;
		}
	}
	// This calculation is wacky, but expressing the shift as a left bit shift
	// within each block lets us use getShiftedBlock.
	Index rightShiftBlocks = (b + N - 1) / N;
	unsigned int leftShiftBits = N * rightShiftBlocks - b;
	// Now (N * rightShiftBlocks - leftShiftBits) == b
	// and 0 <= leftShiftBits < N.
	if (rightShiftBlocks >= a.len + 1) {
		// All of a is guaranteed to be shifted off, even considering the left
		// bit shift.
		len = 0;
		return;
	}
	// Now we're allocating a positive amount.
	// + 1: room for high bits nudged left into another block
	len = a.len + 1 - rightShiftBlocks;
	allocate(len);
	Index i, j;
	for (j = rightShiftBlocks, i = 0; j <= a.len; j++, i++)
		blk[i] = getShiftedBlock(a, j, leftShiftBits);
	// Zap possible leading zero
	if (blk[len - 1] == 0)
		len--;
}

// INCREMENT/DECREMENT OPERATORS

// Prefix increment
void BigUnsigned::operator ++()
{
	Index i;
	bool carry = true;
	for (i = 0; i < len && carry; i++) {
		blk[i]++;
		carry = (blk[i] == 0);
	}
	if (carry) {
		// Allocate and then increase length, as in divideWithRemainder
		allocateAndCopy(len + 1);
		len++;
		blk[i] = 1;
	}
}

// Postfix increment: same as prefix
void BigUnsigned::operator ++(int)
{
	operator ++();
}

// Prefix decrement
void BigUnsigned::operator --()
{
	if (len == 0)
		throw "BigUnsigned::operator --(): Cannot decrement an unsigned zero";
	Index i;
	bool borrow = true;
	for (i = 0; borrow; i++) {
		borrow = (blk[i] == 0);
		blk[i]--;
	}
	// Zap possible leading zero (there can only be one)
	if (blk[len - 1] == 0)
		len--;
}

// Postfix decrement: same as prefix
void BigUnsigned::operator --(int)
{
	operator --();
}
//===================================================================
BigUnsignedInABase::BigUnsignedInABase(const Digit *d, Index l, Base base)
	: NumberlikeArray<Digit>(d, l), base(base) {
	// Check the base
	if (base < 2)
		throw "BigUnsignedInABase::BigUnsignedInABase(const Digit *, Index, Base): The base must be at least 2";

	// Validate the digits.
	for (Index i = 0; i < l; i++)
		if (blk[i] >= base)
			throw "BigUnsignedInABase::BigUnsignedInABase(const Digit *, Index, Base): A digit is too large for the specified base";

	// Eliminate any leading zeros we may have been passed.
	zapLeadingZeros();
}

namespace {
	unsigned int bitLen(unsigned int x) {
		unsigned int len = 0;
		while (x > 0) {
			x >>= 1;
			len++;
		}
		return len;
	}
	unsigned int ceilingDiv(unsigned int a, unsigned int b) {
		return (a + b - 1) / b;
	}
}

BigUnsignedInABase::BigUnsignedInABase(const BigUnsigned &x, Base base) {
	// Check the base
	if (base < 2)
		throw "BigUnsignedInABase(BigUnsigned, Base): The base must be at least 2";
	this->base = base;

	// Get an upper bound on how much space we need
	int maxBitLenOfX = x.getLength() * BigUnsigned::N;
	int minBitsPerDigit = bitLen(base) - 1;
	int maxDigitLenOfX = ceilingDiv(maxBitLenOfX, minBitsPerDigit);
	len = maxDigitLenOfX; // Another change to comply with `staying in bounds'.
	allocate(len); // Get the space

	BigUnsigned x2(x), buBase(base);
	Index digitNum = 0;

	while (!x2.isZero()) {
		// Get last digit.  This is like `lastDigit = x2 % buBase, x2 /= buBase'.
		BigUnsigned lastDigit(x2);
		lastDigit.divideWithRemainder(buBase, x2);
		// Save the digit.
		blk[digitNum] = lastDigit.toUnsignedShort();
		// Move on.  We can't run out of room: we figured it out above.
		digitNum++;
	}

	// Save the actual length.
	len = digitNum;
}

BigUnsignedInABase::operator BigUnsigned() const {
	BigUnsigned ans(0), buBase(base), temp;
	Index digitNum = len;
	while (digitNum > 0) {
		digitNum--;
		temp.multiply(ans, buBase);
		ans.add(temp, BigUnsigned(blk[digitNum]));
	}
	return ans;
}

BigUnsignedInABase::BigUnsignedInABase(const std::string &s, Base base) {
	// Check the base.
	if (base > 36)
		throw "BigUnsignedInABase(std::string, Base): The default string conversion routines use the symbol set 0-9, A-Z and therefore support only up to base 36.  You tried a conversion with a base over 36; write your own string conversion routine.";
	// Save the base.
	// This pattern is seldom seen in C++, but the analogous ``this.'' is common in Java.
	this->base = base;

	// `s.length()' is a `size_t', while `len' is a `NumberlikeArray::Index',
	// also known as an `unsigned int'.  Some compilers warn without this cast.
	len = Index(s.length());
	allocate(len);

	Index digitNum, symbolNumInString;
	for (digitNum = 0; digitNum < len; digitNum++) {
		symbolNumInString = len - 1 - digitNum;
		char theSymbol = s[symbolNumInString];
		if (theSymbol >= '0' && theSymbol <= '9')
			blk[digitNum] = theSymbol - '0';
		else if (theSymbol >= 'A' && theSymbol <= 'Z')
			blk[digitNum] = theSymbol - 'A' + 10;
		else if (theSymbol >= 'a' && theSymbol <= 'z')
			blk[digitNum] = theSymbol - 'a' + 10;
		else
			throw "BigUnsignedInABase(std::string, Base): Bad symbol in input.  Only 0-9, A-Z, a-z are accepted.";

		if (blk[digitNum] >= base)
			throw "BigUnsignedInABase::BigUnsignedInABase(const Digit *, Index, Base): A digit is too large for the specified base";
	}
	zapLeadingZeros();
}

BigUnsignedInABase::operator std::string() const {
	if (base > 36)
		throw "BigUnsignedInABase ==> std::string: The default string conversion routines use the symbol set 0-9, A-Z and therefore support only up to base 36.  You tried a conversion with a base over 36; write your own string conversion routine.";
	if (len == 0)
		return std::string("0");
	// Some compilers don't have push_back, so use a char * buffer instead.
	char *s = new char[len + 1];
	s[len] = '\0';
	Index digitNum, symbolNumInString;
	for (symbolNumInString = 0; symbolNumInString < len; symbolNumInString++) {
		digitNum = len - 1 - symbolNumInString;
		Digit theDigit = blk[digitNum];
		if (theDigit < 10)
			s[symbolNumInString] = char('0' + theDigit);
		else
			s[symbolNumInString] = char('A' + theDigit - 10);
	}
	std::string s2(s);
	delete [] s;
	return s2;
}




void BigInteger::operator =(const BigInteger &x)
{
	// Calls like a = a have no effect
	if (this == &x)
		return;
	// Copy sign
	sign = x.sign;
	// Copy the rest
	mag = x.mag;
}

BigInteger::BigInteger(const Blk *b, Index blen, Sign s) :
		mag(b, blen)
{
	switch (s) {
	case zero:
		if (!mag.isZero())
			throw "BigInteger::BigInteger(const Blk *, Index, Sign): Cannot use a sign of zero with a nonzero magnitude";
		sign = zero;
		break;
	case positive:
	case negative:
		// If the magnitude is zero, force the sign to zero.
		sign = mag.isZero() ? zero : s;
		break;
	default:
		/* g++ seems to be optimizing out this case on the assumption
		 * that the sign is a valid member of the enumeration.  Oh well. */
		throw "BigInteger::BigInteger(const Blk *, Index, Sign): Invalid sign";
	}
}

BigInteger::BigInteger(const BigUnsigned &x, Sign s) :
		mag(x)
{
	switch (s) {
	case zero:
		if (!mag.isZero())
			throw "BigInteger::BigInteger(const BigUnsigned &, Sign): Cannot use a sign of zero with a nonzero magnitude";
		sign = zero;
		break;
	case positive:
	case negative:
		// If the magnitude is zero, force the sign to zero.
		sign = mag.isZero() ? zero : s;
		break;
	default:
		/* g++ seems to be optimizing out this case on the assumption
		 * that the sign is a valid member of the enumeration.  Oh well. */
		throw "BigInteger::BigInteger(const BigUnsigned &, Sign): Invalid sign";
	}
}

/* CONSTRUCTION FROM PRIMITIVE INTEGERS
 * Same idea as in BigUnsigned.cc, except that negative input results in a
 * negative BigInteger instead of an exception. */

// Done longhand to let us use initialization.
BigInteger::BigInteger(unsigned long x) :
		mag(x)
{
	sign = mag.isZero() ? zero : positive;
}
BigInteger::BigInteger(unsigned int x) :
		mag(x)
{
	sign = mag.isZero() ? zero : positive;
}
BigInteger::BigInteger(unsigned short x) :
		mag(x)
{
	sign = mag.isZero() ? zero : positive;
}

// For signed input, determine the desired magnitude and sign separately.

namespace
{
template<class X, class UX>
BigInteger::Blk magOf(X x)
{
	/* UX(...) cast needed to stop short(-2^15), which negates to
	 * itself, from sign-extending in the conversion to Blk. */
	return BigInteger::Blk(x < 0 ? UX(-x) : x);
}
template<class X>
BigInteger::Sign signOf(X x)
{
	return (x == 0) ? BigInteger::zero :
			(x > 0) ? BigInteger::positive : BigInteger::negative;
}
}

BigInteger::BigInteger(long x) :
		sign(signOf(x)), mag(magOf<long, unsigned long>(x))
{
}
BigInteger::BigInteger(int x) :
		sign(signOf(x)), mag(magOf<int, unsigned int>(x))
{
}
BigInteger::BigInteger(short x) :
		sign(signOf(x)), mag(magOf<short, unsigned short>(x))
{
}

// CONVERSION TO PRIMITIVE INTEGERS

/* Reuse BigUnsigned's conversion to an unsigned primitive integer.
 * The friend is a separate function rather than
 * BigInteger::convertToUnsignedPrimitive to avoid requiring BigUnsigned to
 * declare BigInteger. */
template<class X>
inline X convertBigUnsignedToPrimitiveAccess(const BigUnsigned &a)
{
	return a.convertToPrimitive<X>();
}

template<class X>
X BigInteger::convertToUnsignedPrimitive() const
{
	if (sign == negative)
		throw "BigInteger::to<Primitive>: "
				"Cannot convert a negative integer to an unsigned type";
	else
		return convertBigUnsignedToPrimitiveAccess<X>(mag);
}

/* Similar to BigUnsigned::convertToPrimitive, but split into two cases for
 * nonnegative and negative numbers. */
template<class X, class UX>
X BigInteger::convertToSignedPrimitive() const
{
	if (sign == zero)
		return 0;
	else if (mag.getLength() == 1) {
		// The single block might fit in an X.  Try the conversion.
		Blk b = mag.getBlock(0);
		if (sign == positive) {
			X x = X(b);
			if (x >= 0 && Blk(x) == b)
				return x;
		} else {
			X x = -X(b);
			/* UX(...) needed to avoid rejecting conversion of
			 * -2^15 to a short. */
			if (x < 0 && Blk(UX(-x)) == b)
				return x;
		}
		// Otherwise fall through.
	}
	throw "BigInteger::to<Primitive>: "
			"Value is too big to fit in the requested type";
}

unsigned long BigInteger::toUnsignedLong() const
{
	return convertToUnsignedPrimitive<unsigned long>();
}
unsigned int BigInteger::toUnsignedInt() const
{
	return convertToUnsignedPrimitive<unsigned int>();
}
unsigned short BigInteger::toUnsignedShort() const
{
	return convertToUnsignedPrimitive<unsigned short>();
}
long BigInteger::toLong() const
{
	return convertToSignedPrimitive<long, unsigned long>();
}
int BigInteger::toInt() const
{
	return convertToSignedPrimitive<int, unsigned int>();
}
short BigInteger::toShort() const
{
	return convertToSignedPrimitive<short, unsigned short>();
}

// COMPARISON
BigInteger::CmpRes BigInteger::compareTo(const BigInteger &x) const
{
	// A greater sign implies a greater number
	if (sign < x.sign)
		return less;
	else if (sign > x.sign)
		return greater;
	else
		switch (sign) {
		// If the signs are the same...
		case zero:
			return equal; // Two zeros are equal
		case positive:
			// Compare the magnitudes
			return mag.compareTo(x.mag);
		case negative:
			// Compare the magnitudes, but return the opposite result
			return CmpRes(-mag.compareTo(x.mag));
		default:
			throw "BigInteger internal error";
		}
}

/* COPY-LESS OPERATIONS
 * These do some messing around to determine the sign of the result,
 * then call one of BigUnsigned's copy-less operations. */
// See remarks about aliased calls in BigUnsigned.cc .
#define DTRT_ALIASED_BI(cond, op) \
	if (cond) { \
		BigInteger tmpThis; \
		tmpThis.op; \
		*this = tmpThis; \
		return; \
	}
void BigInteger::add(const BigInteger &a, const BigInteger &b)
{
	DTRT_ALIASED_BI(this == &a || this == &b, add(a, b));
	// If one argument is zero, copy the other.
	if (a.sign == zero)
		operator =(b);
	else if (b.sign == zero)
		operator =(a);
	// If the arguments have the same sign, take the
	// common sign and add their magnitudes.
	else if (a.sign == b.sign) {
		sign = a.sign;
		mag.add(a.mag, b.mag);
	} else {
		// Otherwise, their magnitudes must be compared.
		switch (a.mag.compareTo(b.mag)) {
		case equal:
			// If their magnitudes are the same, copy zero.
			mag = 0;
			sign = zero;
			break;
			// Otherwise, take the sign of the greater, and subtract
			// the lesser magnitude from the greater magnitude.
		case greater:
			sign = a.sign;
			mag.subtract(a.mag, b.mag);
			break;
		case less:
			sign = b.sign;
			mag.subtract(b.mag, a.mag);
			break;
		}
	}
}

void BigInteger::subtract(const BigInteger &a, const BigInteger &b)
{
	// Notice that this routine is identical to BigInteger::add,
	// if one replaces b.sign by its opposite.
	DTRT_ALIASED_BI(this == &a || this == &b, subtract(a, b));
	// If a is zero, copy b and flip its sign.  If b is zero, copy a.
	if (a.sign == zero) {
		mag = b.mag;
		// Take the negative of _b_'s, sign, not ours.
		// Bug pointed out by Sam Larkin on 2005.03.30.
		sign = Sign(-b.sign);
	} else if (b.sign == zero)
		operator =(a);
	// If their signs differ, take a.sign and add the magnitudes.
	else if (a.sign != b.sign) {
		sign = a.sign;
		mag.add(a.mag, b.mag);
	} else {
		// Otherwise, their magnitudes must be compared.
		switch (a.mag.compareTo(b.mag)) {
		// If their magnitudes are the same, copy zero.
		case equal:
			mag = 0;
			sign = zero;
			break;
			// If a's magnitude is greater, take a.sign and
			// subtract a from b.
		case greater:
			sign = a.sign;
			mag.subtract(a.mag, b.mag);
			break;
			// If b's magnitude is greater, take the opposite
			// of b.sign and subtract b from a.
		case less:
			sign = Sign(-b.sign);
			mag.subtract(b.mag, a.mag);
			break;
		}
	}
}

void BigInteger::multiply(const BigInteger &a, const BigInteger &b)
{
	DTRT_ALIASED_BI(this == &a || this == &b, multiply(a, b));
	// If one object is zero, copy zero and return.
	if (a.sign == zero || b.sign == zero) {
		sign = zero;
		mag = 0;
		return;
	}
	// If the signs of the arguments are the same, the result
	// is positive, otherwise it is negative.
	sign = (a.sign == b.sign) ? positive : negative;
	// Multiply the magnitudes.
	mag.multiply(a.mag, b.mag);
}

/*
 * DIVISION WITH REMAINDER
 * Please read the comments before the definition of
 * `BigUnsigned::divideWithRemainder' in `BigUnsigned.cc' for lots of
 * information you should know before reading this function.
 *
 * Following Knuth, I decree that x / y is to be
 * 0 if y==0 and floor(real-number x / y) if y!=0.
 * Then x % y shall be x - y*(integer x / y).
 *
 * Note that x = y * (x / y) + (x % y) always holds.
 * In addition, (x % y) is from 0 to y - 1 if y > 0,
 * and from -(|y| - 1) to 0 if y < 0.  (x % y) = x if y = 0.
 *
 * Examples: (q = a / b, r = a % b)
 *	a	b	q	r
 *	===	===	===	===
 *	4	3	1	1
 *	-4	3	-2	2
 *	4	-3	-2	-2
 *	-4	-3	1	-1
 */
void BigInteger::divideWithRemainder(const BigInteger &b, BigInteger &q)
{
	// Defend against aliased calls;
	// same idea as in BigUnsigned::divideWithRemainder .
	if (this == &q)
		throw "BigInteger::divideWithRemainder: Cannot write quotient and remainder into the same variable";
	if (this == &b || &q == &b) {
		BigInteger tmpB(b);
		divideWithRemainder(tmpB, q);
		return;
	}

	// Division by zero gives quotient 0 and remainder *this
	if (b.sign == zero) {
		q.mag = 0;
		q.sign = zero;
		return;
	}
	// 0 / b gives quotient 0 and remainder 0
	if (sign == zero) {
		q.mag = 0;
		q.sign = zero;
		return;
	}

	// Here *this != 0, b != 0.

	// Do the operands have the same sign?
	if (sign == b.sign) {
		// Yes: easy case.  Quotient is zero or positive.
		q.sign = positive;
	} else {
		// No: harder case. Quotient is negative.
		q.sign = negative;
		// Decrease the magnitude of the dividend by one.
		mag--;
		/*
		 * We tinker with the dividend before and with the
		 * quotient and remainder after so that the result
		 * comes out right.  To see why it works, consider the following
		 * list of examples, where A is the magnitude-decreased
		 * a, Q and R are the results of BigUnsigned division
		 * with remainder on A and |b|, and q and r are the
		 * final results we want:
		 *
		 *	a	A	b	Q	R	q	r
		 *	-3	-2	3	0	2	-1	0
		 *	-4	-3	3	1	0	-2	2
		 *	-5	-4	3	1	1	-2	1
		 *	-6	-5	3	1	2	-2	0
		 *
		 * It appears that we need a total of 3 corrections:
		 * Decrease the magnitude of a to get A.  Increase the
		 * magnitude of Q to get q (and make it negative).
		 * Find r = (b - 1) - R and give it the desired sign.
		 */
	}

	// Divide the magnitudes.
	mag.divideWithRemainder(b.mag, q.mag);

	if (sign != b.sign) {
		// More for the harder case (as described):
		// Increase the magnitude of the quotient by one.
		q.mag++;
		// Modify the remainder.
		mag.subtract(b.mag, mag);
		mag--;
	}

	// Sign of the remainder is always the sign of the divisor b.
	sign = b.sign;

	// Set signs to zero as necessary.  (Thanks David Allen!)
	if (mag.isZero())
		sign = zero;
	if (q.mag.isZero())
		q.sign = zero;

	// WHEW!!!
}

// Negation
void BigInteger::negate(const BigInteger &a)
{
	DTRT_ALIASED_BI(this == &a, negate(a));
	// Copy a's magnitude
	mag = a.mag;
	// Copy the opposite of a.sign
	sign = Sign(-a.sign);
}

// INCREMENT/DECREMENT OPERATORS

// Prefix increment
void BigInteger::operator ++()
{
	if (sign == negative) {
		mag--;
		if (mag == 0)
			sign = zero;
	} else {
		mag++;
		sign = positive; // if not already
	}
}

// Postfix increment: same as prefix
void BigInteger::operator ++(int)
{
	operator ++();
}

// Prefix decrement
void BigInteger::operator --()
{
	if (sign == positive) {
		mag--;
		if (mag == 0)
			sign = zero;
	} else {
		mag++;
		sign = negative;
	}
}

// Postfix decrement: same as prefix
void BigInteger::operator --(int)
{
	operator --();
}

std::string bigUnsignedToString(const BigUnsigned &x) {
	return std::string(BigUnsignedInABase(x, 10));
}

std::string bigIntegerToString(const BigInteger &x) {
	return (x.getSign() == BigInteger::negative)
		? (std::string("-") + bigUnsignedToString(x.getMagnitude()))
		: (bigUnsignedToString(x.getMagnitude()));
}

BigUnsigned stringToBigUnsigned(const std::string &s) {
	return BigUnsigned(BigUnsignedInABase(s, 10));
}

BigInteger stringToBigInteger(const std::string &s) {
	// Recognize a sign followed by a BigUnsigned.
	return (s[0] == '-') ? BigInteger(stringToBigUnsigned(s.substr(1, s.length() - 1)), BigInteger::negative)
		: (s[0] == '+') ? BigInteger(stringToBigUnsigned(s.substr(1, s.length() - 1)))
		: BigInteger(stringToBigUnsigned(s));
}

std::ostream &operator <<(std::ostream &os, const BigUnsigned &x) {
	BigUnsignedInABase::Base base;
	long osFlags = os.flags();
	if (osFlags & os.dec)
		base = 10;
	else if (osFlags & os.hex) {
		base = 16;
		if (osFlags & os.showbase)
			os << "0x";
	} else if (osFlags & os.oct) {
		base = 8;
		if (osFlags & os.showbase)
			os << '0';
	} else
		throw "std::ostream << BigUnsigned: Could not determine the desired base from output-stream flags";
	std::string s = std::string(BigUnsignedInABase(x, base));
	os << s;
	return os;
}

std::ostream &operator <<(std::ostream &os, const BigInteger &x) {
	if (x.getSign() == BigInteger::negative)
		os << '-';
	os << x.getMagnitude();
	return os;
}

//BigUnsigned gcd(BigUnsigned a, BigUnsigned b) {
//	BigUnsigned trash;
	// Neat in-place alternating technique.
//	for (;;) {
//		if (b.isZero())
//			return a;
//		a.divideWithRemainder(b, trash);
//		if (a.isZero())
//			return b;
//		b.divideWithRemainder(a, trash);
//	}
//}

void extendedEuclidean(BigInteger m, BigInteger n,
		BigInteger &g, BigInteger &r, BigInteger &s) {
	if (&g == &r || &g == &s || &r == &s)
		throw "BigInteger extendedEuclidean: Outputs are aliased";
	BigInteger r1(1), s1(0), r2(0), s2(1), q;
	/* Invariants:
	 * r1*m(orig) + s1*n(orig) == m(current)
	 * r2*m(orig) + s2*n(orig) == n(current) */
	for (;;) {
		if (n.isZero()) {
			r = r1; s = s1; g = m;
			return;
		}
		// Subtract q times the second invariant from the first invariant.
		m.divideWithRemainder(n, q);
		r1 -= q*r2; s1 -= q*s2;

		if (m.isZero()) {
			r = r2; s = s2; g = n;
			return;
		}
		// Subtract q times the first invariant from the second invariant.
		n.divideWithRemainder(m, q);
		r2 -= q*r1; s2 -= q*s1;
	}
}

BigUnsigned modinv(const BigInteger &x, const BigUnsigned &n) {
	BigInteger g, r, s;
	extendedEuclidean(x, n, g, r, s);
	if (g == 1)
		// r*x + s*n == 1, so r*x === 1 (mod n), so r is the answer.
		return (r % n).getMagnitude(); // (r % n) will be nonnegative
	else
		throw "BigInteger modinv: x and n have a common factor";
}

BigUnsigned modexp(const BigInteger &base, const BigUnsigned &exponent,
		const BigUnsigned &modulus) {
	BigUnsigned ans = 1, base2 = (base % modulus).getMagnitude();
	BigUnsigned::Index i = exponent.bitLength();
	// For each bit of the exponent, most to least significant...
	while (i > 0) {
		i--;
		// Square.
		ans *= ans;
		ans %= modulus;
		// And multiply if the bit is a 1.
		if (exponent.getBit(i)) {
			ans *= base2;
			ans %= modulus;
		}
	}
	return ans;
}


}

