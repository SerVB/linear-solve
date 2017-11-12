class Matrix:
    def __init__(self, source=[[]]):
        self.storage = []
        for row in source:
            nowRow = []
            for cell in row:
                nowRow.append(cell)
            self.storage.append(nowRow)

    def __str__(self):
        rowStrrings = []
        for row in self.storage:
            rowStrrings.append("\t".join(map(str, row)))
        return "\n".join(rowStrrings)

    def size(self):
        # Row count, column count
        return (len(self.storage), len(self.storage[0]))

    def __mul__(self, object):
        if isinstance(object, Matrix):
            return self.mulToMatrix(object)
        else:
            return self.mulToScalar(object)

    def mulToScalar(self, number):
        newStorage = []
        for row in self.storage:
            newRow = []
            for cell in row:
                newRow.append(cell * number)
            newStorage.append(newRow)
        return Matrix(newStorage)

    def getColumn(self, number):
        for row in self.storage:
            yield row[number]

    def mulToMatrix(self, matr):
        if self.size()[1] != matr.size()[0]:
            raise MatrixError(self, matr)
        newStorage = []
        for row in self.storage:
            newRow = []
            for columnIdx in range(matr.size()[1]):
                summ = sum(
                    map(
                        lambda x: x[0] * x[1],
                        zip(row, matr.getColumn(columnIdx))
                    )
                )
                newRow.append(summ)
            newStorage.append(newRow)
        return Matrix(newStorage)

    __rmul__ = __mul__

    def __add__(self, other):
        if self.size() != other.size():
            raise MatrixError(self, other)

        newStorage = []
        for row in self.storage:
            newRow = []
            for cell in row:
                newRow.append(cell)
            newStorage.append(newRow)

        for i in range(len(newStorage)):
            for j in range(len(newStorage[0])):
                newStorage[i][j] += other.storage[i][j]

        return Matrix(newStorage)

    def transpose(self):
        newStorage = []
        for columnIdx in range(len(self.storage[0])):
            newRow = list(self.getColumn(columnIdx))
            newStorage.append(newRow)
        self.storage = newStorage
        return self

    @staticmethod
    def transposed(matr):
        newMatr = Matrix(matr.storage)
        return newMatr.transpose()

    def getDeterminant(self):
        # Should be in SquareMatrix class!
        assert self.size()[0] == self.size()[1]
        ans = 0
        #for perm in permutations(range(self.size()[0])):
        for perm in permutations(self.size()[0]):
            now = 1
            for i in range(self.size()[0]):
                now *= self.storage[i][perm[i]]
            ans += sgn(perm) * now
        return ans

    def solve(self, vect):
        return self.gauss(vect)
#         return kramer(vect)

    def gauss(self, vect):
        newStorage = []
        for row in self.storage:
            newRow = []
            for cell in row:
                newRow.append(cell)
            newStorage.append(newRow)
        for i in range(len(vect)):
            newStorage[i].append(vect[i])

        for columnIdx in range(len(newStorage[0]) - 1):
            newStorage.sort(reverse=True)
            if newStorage[columnIdx][columnIdx] == 0:
                assert False

            div = newStorage[columnIdx][columnIdx]

            for i in range(columnIdx, len(newStorage[0])):
                newStorage[columnIdx][i] /= div

            for j in range(0, len(newStorage)):
                if j != columnIdx:
                    mul = newStorage[j][columnIdx] / \
                        newStorage[columnIdx][columnIdx]
                    for i in range(len(newStorage[0])):
                        newStorage[j][i] -= mul * newStorage[columnIdx][i]

        ans = Matrix(newStorage)
        return list(ans.getColumn(-1))

    def kramer(self, vect):
        ans = []
        det = self.getDeterminant()
        nowMatr = Matrix(self.storage)

        for i in range(len(vect)):
            for rowIdx in range(len(vect)):
                nowMatr.storage[rowIdx][i] = vect[i]
            nowDet = nowMatr.getDeterminant()

            for rowIdx in range(len(vect)):
                nowMatr.storage[rowIdx][i] = self.storage[rowIdx][i]

            ans.append(nowDet / det)

        return ans


class SquareMatrix(Matrix):

    @staticmethod
    def getId(n):
        # Returns diag(1, 1, ..., 1)
        newStorage = []
        for i in range(n):
            newRow = []
            for j in range(n):
                if i == j:
                    newRow.append(1)
                else:
                    newRow.append(0)
            newStorage.append(newRow)
        return SquareMatrix(newStorage)

    def __pow__(self, n):
        if n < 0:
            assert False
        if n == 0:
            return SquareMatrix.getId(self.size()[0])

        if n % 2 == 0:
            half = SquareMatrix(self.storage) ** (n // 2)
            return half * half

        almostHalf = SquareMatrix(self.storage) ** ((n - 1) // 2)
        return almostHalf * almostHalf * SquareMatrix(self.storage)


class MatrixError(BaseException):
    def __init__(self, matrix1, matrix2):
        self.matrix1 = matrix1
        self.matrix2 = matrix2


class Poly:
    #  coef means dict "degree: value"
    def __init__(self, coefs={0: 0}):
        self.coefs = coefs
        maxDeg = max(coefs.keys())
        for deg in range(maxDeg):
            if deg not in self.coefs:
                self.coefs[deg] = 0

    def __str__(self):
        ans = ""
        for deg, val in self.coefs.items():
            ans += str(val) + "x" + str(deg) + " "
        return ans.rstrip()

    def __add__(self, other):
        if not isinstance(other, Poly):
            other = Poly({0: other})

        if len(other.coefs) > len(self.coefs):
            large = dict(other.coefs)
            small = dict(self.coefs)
        else:
            large = dict(self.coefs)
            small = dict(other.coefs)

        for deg, val in large.items():
            large[deg] += small.get(deg, 0)
        return Poly(large)

    __radd__ = __add__

    def __mul__(self, other):
        if not isinstance(other, Poly):
            other = Poly({0: other})

        newCoefs = dict()

        for deg1, val1 in self.coefs.items():
            for deg2, val2 in other.coefs.items():
                newCoefs[deg1 + deg2] = newCoefs.get(deg1 + deg2, 0) + \
                    val1 * val2

        return Poly(newCoefs)

    __rmul__ = __mul__


class Substitution:
    #  map means dict "i: s[i]" 1..n
    def __init__(self, map):
        self.map = map

    def __mul__(self, other):
        newMap = dict()
        for i in self.map.keys():
            newMap[i] = self.map[other.map[i]]
        return Substitution(newMap)

    def __str__(self):
        ans = list()
        for k, v in self.map.items():
            ans.append(str(k) + ": " + str(v))

        return ", ".join(ans)

    def __pow__(self, n):
        if n == -1:
            newMap = dict()
            for k, v in self.map.items():
                newMap[v] = k
            return Substitution(newMap)
        if n == 0:
            newMap = dict()
            for k, v in self.map.items():
                newMap[k] = k
            return Substitution(newMap)

        initial = Substitution(self.map)
        current = Substitution(self.map)

        while n > 1:
            current *= initial
            n -= 1

        return current


def sgn(permut):
    num = 0
    for i in range(len(permut)):
        for j in range(i + 1, len(permut)):
            if permut[j] < permut[i]:
                num += 1
    return 1 if num % 2 == 0 else -1


def nextPerm(listPerm):
    # Find longest non-increasing suffix
    i = len(listPerm) - 1
    while i > 0 and listPerm[i - 1] >= listPerm[i]:
        i -= 1
    # Now i is the head index of the suffix

    # Are we at the last permutation already?
    if i <= 0:
        return

    # Let listPerm[i - 1] be the pivot
    # Find rightmost element that exceeds the pivot
    j = len(listPerm) - 1
    while listPerm[j] <= listPerm[i - 1]:
        j -= 1
    # Now the value listPerm[j] will become the new pivot
    assert j >= i

    # Swap the pivot with j
    temp = listPerm[i - 1]
    listPerm[i - 1] = listPerm[j]
    listPerm[j] = temp

    # Reverse the suffix
    j = len(listPerm) - 1
    while i < j:
        temp = listPerm[i]
        listPerm[i] = listPerm[j]
        listPerm[j] = temp
        i += 1
        j -= 1


def permutations(n):
    target = tuple(range(n)[::-1])
    now = list(range(n))
    while tuple(now) != target:
        yield tuple(now)
        nextPerm(now)
    yield target


#  Example Task 1
print("Task 1:")  # ans is 1: 8, 2: 6, 3: 5, 4: 2, 5: 7, 6: 3, 7: 1, 8: 4
a = Substitution({1: 5, 2: 8, 3: 4, 4: 3, 5: 2, 6: 1, 7: 6, 8: 7}) ** 19
b = Substitution({1: 4, 2: 5, 3: 1, 4: 3, 5: 6, 6: 8, 7: 7, 8: 2}) ** -1
left = (a * b) ** 159
right = Substitution({1: 5, 2: 6, 3: 4, 4: 2, 5: 7, 6: 3, 7: 8, 8: 1})

perm = [1, 2, 3, 4, 5, 6, 7, 8]
target = tuple(perm[::-1])
while tuple(perm) != target:
    newMap = dict()
    for i in range(8):
        newMap[i + 1] = perm[i]
    x = Substitution(newMap)
    res = left * x
    if res.map == right.map:
        print(x)
    nextPerm(perm)
#  Check the last one:
newMap = dict()
for i in range(8):
    newMap[i + 1] = target[i]
x = Substitution(newMap)
res = left * x
if res.map == right.map:
    print(x)
#  Example Task 1

#  Example Task 2
perm = []
for i in range(1, 126 + 1):
    perm.append(i + 21)
for i in range(127, 147 + 1):
    perm.append(i - 126)
i = 1
print("Task 2 sequence:")
for elem in perm:
    print("%d:%d" % (i, elem), end=" ")
    i += 1

print("\nTask 2:", sgn(perm))  # 1
#  Example Task 2

#  Example Task 3
l = [
    [0, 0, 3, 1, 0, 0],
    [0, 2, 2, 2, 3, 0],
    [0, 1, 3, 0, 4, 3],
    [0, 3, 1, 4, 3, 0],
    [1, 0, 0, 6, 3, 4],
    [1, 7, 0, 2, 2, 3]
    ]
m = Matrix(l)

print("Task 3:", m.getDeterminant())  # 569
#  Example Task 3

#  Example Task 4
x = Poly({1: 1})
l = [
    [7, 8, 5, 3, 9, x, 1],
    [8, 10, 7, 5, x, 4, 3],
    [5, 7, 5, x, 9, 2, 1],
    [3, 5, x, 5, 0, 4, 3],
    [9, x, 9, 0, 9, 2, 1],
    [x, 4, 2, 4, 2, 4, x],
    [1, 3, 1, 3, 1, x, 1]
    ]
m = Matrix(l)
# -6x6 80x5 229x4 -8132x3 44862x2 -97784x1 77064x0
print("Task 4:", m.getDeterminant())
#  Example Task 4
