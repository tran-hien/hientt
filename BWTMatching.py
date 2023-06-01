class SuffixArray:
    '''
    Tạo mảng hậu tố bằng thuật toán SA-IS
    '''
    S_type = ord("S")
    L_type = ord("L")
    def __init__(self, text, alphabetSize=256):
        self.suffix_array = self.makeSuffixArray(text, alphabetSize)

    def catergorize(self, Text):
        t = bytearray(len(Text) + 1)  #type array
        t[-1] = self.S_type
        t[-2] = self.L_type
        #Duyệt Text từ phải sang trái
        for i in range(len(Text)-2,-1,-1):
            if Text[i] < Text[i+1] or (Text[i] == Text[i+1] and t[i+1] == self.S_type):
                t[i] = self.S_type
            else:
                t[i] = self.L_type         
        return t 

    #Hàm check ký tự thứ index trong Text có phải loại LSM ko
    def isLMSChar(self, index, t):
        if index == 0:
            return False
        if t[index] == self.S_type and t[index-1] == self.L_type:
            return True
        return False

    
    def lmsSubstringsAreEqual(self, string, typemap, offsetA, offsetB):
        """
        Return True nếu LMS substrings tại offsetA và offsetB bằng nhau.
        """
        if offsetA == len(string) or offsetB == len(string):
            return False

        i = 0
        while True:
            aIsLMS = self.isLMSChar(i + offsetA, typemap)
            bIsLMS = self.isLMSChar(i + offsetB, typemap)

        # Nếu tìm được vị trí bắt đầu của LMS substrings tiếp theo
            if (i > 0 and aIsLMS and bIsLMS):
                return True

            if aIsLMS != bIsLMS:
                return False

            if string[i + offsetA] != string[i + offsetB]:
                return False

            i += 1


    def findBucketSize(self, string, alphabetSize=256):
        B = [0] * alphabetSize
        for char in string:
            B[char] += 1  
        return B

    def findBucketHeads(self, bucketSize):
        start = 1
        res = []
        for size in bucketSize:
            res.append(start)
            start += size        
        return res

    def findBucketTails(self, bucketSize):
        end = 1
        res = []
        for size in bucketSize:
            end += size
            res.append(end - 1)  
        return res


    def makeSuffixArrayByInducedSorting(self, string, alphabetSize):
        """
        Hàm chính: Tìm mảng hậu tố bằng thuật toán SA-IS
        """

        # Tìm loại S/L của các hậu tố
        typemap = self.catergorize(string)

        #Tìm cỡ của các bucket
        bucketSizes = self.findBucketSize(string, alphabetSize)

        #Tạo mảng hậu tố tạm, insert vị trí các hậu tố loại LMS
        guessedSuffixArray = self.guessLMSSort(string, bucketSizes, typemap)
        #induced-Sort các hậu tố loại L, S
        self.inducedSortL(string, guessedSuffixArray, bucketSizes, typemap)
        self.inducedSortS(string, guessedSuffixArray, bucketSizes, typemap)

        #Tạo T1 mới, bảng chữ cái mới 
        summaryString, summaryAlphabetSize, summarySuffixOffsets = self.summarySA(string, guessedSuffixArray, typemap)

        # Gọi đến hàm đệ quy tìm T1 và bảng chữ cái mới để tìm được mảng hậu tố rút gọn
        summarySuffixArray = self.recusivelyMakeSA(
            summaryString,
            summaryAlphabetSize,
        )

        # xếp vị trí chính xác của các hậu tố LMS vào mảng hậu tố cuối cùng
        result = self.accurateLMSSort(string, bucketSizes, typemap,
                summarySuffixArray, summarySuffixOffsets)

        # Xếp vị trí các hậu tố loại L và loại S vào mảng hậu tố cuối cùng
        self.inducedSortL(string, result, bucketSizes, typemap)
        self.inducedSortS(string, result, bucketSizes, typemap)

        return result

    def guessLMSSort(self, string, bucketSizes, typemap):
        guessedSuffixArray = [-1] * (len(string) + 1)
        
        bucketTails = self.findBucketTails(bucketSizes)
        
        for i in range(len(string)):
            if not self.isLMSChar(i, typemap):
                #khong phai diem bat dau cua hau to LMS
                continue
            
            #string[i] la thuoc bucket nao
            bucketIndex = string[i] #=c/a/b/b/ sau khi da duoc encode, tuc la no la kieu integer
            
            #thêm biến trên vào mảng tạm SA
            guessedSuffixArray[bucketTails[bucketIndex]] = i
            #vd: guessSuffixArray[bucketTails[0(~a)]] = guessSuffixArray[2]
            #Tức đưa ký tự bắt đầu của chuỗi hậu tố LMS vào cuối của bucket tương ứng. 
            # Ký tự đầu là gì thì nó thuộc bucket đấy
            
            # Dịch đuôi bucket sang trái 1 đơn vị
            bucketTails[bucketIndex] -= 1
        
        #bucket đầu là của ký tự $, trong code này là của ký tự cuối của String là rỗng thôi  
        guessedSuffixArray[0] = len(string)
        
        #showSuffixArray(guessedSuffixArray)
        
        return guessedSuffixArray

#--------------------Induced Sorting L-type--------------------------
    def inducedSortL(self, text, guessedSA, bucketSizes, typemap):
        bucketHeads = self.findBucketHeads(bucketSizes)
        
        for i in range(len(guessedSA)):
            if guessedSA[i] == -1:
                continue
            
            
            if typemap[guessedSA[i]-1] != self.L_type:
                continue
            
            bucketIndex = text[guessedSA[i]-1]
            guessedSA[bucketHeads[bucketIndex]] = guessedSA[i]-1
                
            bucketHeads[bucketIndex] += 1

            #showSuffixArray(guessedSA, i)
        return guessedSA


#----------------Induced Sorting S----------------------
    def inducedSortS(self, text, guessedSA, bucketSizes, typemap):
        bucketTails = self.findBucketTails(bucketSizes)
        
        for i in range(len(guessedSA)-1,-1,-1):
            j = guessedSA[i]-1
            
            if j < 0:
                continue
            
            if typemap[j] == self.S_type:
                bucketIndex = text[j]
                
                guessedSA[bucketTails[bucketIndex]] = j
                
                bucketTails[bucketIndex] -=1
                
                #showSuffixArray(guessedSA, i)
            
        return guessedSA
 
    #Hàm trả về chuỗi mới(T1) và Bảng chữ cái mới        
    def summarySA(self, text, guessedSA, typemap):
        #Mảng lưu trữ tên mới của các chuỗi con LMS
        lmsNames = [-1] * (len(text) + 1)
        
        #Tên hiện tại của LMS-substring
        currentName = 0
        
        #Ta biết các LMS-substring sau khi sắp xếp thì $ luôn đứng đầu nên ta đặt tên nó là 0
        lmsNames[guessedSA[0]] = currentName
        lastLMSSubstringOffset = guessedSA[0]
        
        #Duyệt mảng SA:
        for i in range(1, len(guessedSA)):
            #Bỏ qua nếu ký tự thứ guessedSA[i] của text không phải là ký tự bắt đầu của chuỗi LMS
            if not self.isLMSChar(guessedSA[i], typemap):
                continue
            
            #Tăng tên lên 1 đơn vị nếu LMS tiếp theo không giống cái trước = > không cùng tên
            if not self.lmsSubstringsAreEqual(text, typemap, lastLMSSubstringOffset, guessedSA[i]):
                currentName += 1
            
            #Đặt lại biến này vì LMS đã thay đổi
            lastLMSSubstringOffset = guessedSA[i]
            
            #Gán tên mới của LMS vào mảng tên với vị trí tương ứng vị trí của LMS trong text
            lmsNames[lastLMSSubstringOffset] = currentName
            
            summarSuffixOffsets = []
            T1 = []
    
        for i in range(len(lmsNames)):
            if lmsNames[i] != -1:
                summarSuffixOffsets.append(i)
                T1.append(lmsNames[i])
                    
        newAlphabetSize = currentName + 1
        
        return T1, newAlphabetSize, summarSuffixOffsets     
        

    #-----------Bước 5: Đệ quy---------------------
    def recusivelyMakeSA(self, T1, newAlphabetSize):
        #Nếu chuỗi mới T1 gồm các ký tự không lặp lại => base case
        if newAlphabetSize == len(T1):
            SAfinal = [-1] * (len(T1) + 1)
            
            #Phần tử đầu của SA luôn thuộc về $
            SAfinal[0] = len(T1)
            
            for i in range(len(T1)):
                y= T1[i]
                SAfinal[y+1] = i
                
        else:
            SAfinal = self.makeSuffixArrayByInducedSorting(T1, newAlphabetSize)
        return SAfinal

        
# def showSuffixArray(arr, pos=None):
#     print(" ".join("%02d" % each for each in arr))

#     if pos is not None:
#         print(" ".join(
#         "^^" if each == pos else "  "
#                 for each in range(len(arr))
#         ))

    def accurateLMSSort(self, string, bucketSizes, typemap,
        summarySuffixArray, summarySuffixOffsets):
        """
        Tạo mảng hậu tố chứa các vị trí các hậu tố LMS đúng chỗ
        """
        suffixOffsets = [-1] * (len(string) + 1)

        # Trước đây ta thêm các hậu tố vào cuối mỗi bucket tương ứng, vì thế để 
        # giữ chúng ở vị trú đúng ta sẽ phải duyệt summarySuffixArray ngược
        bucketTails = self.findBucketTails(bucketSizes)
        for i in range(len(summarySuffixArray)-1, 1, -1):
            stringIndex = summarySuffixOffsets[summarySuffixArray[i]]

            # Tìm bucket
            bucketIndex = string[stringIndex]
            # Thêm hậu tố vào cuối bucket
            suffixOffsets[bucketTails[bucketIndex]] = stringIndex
            # Lùi đuôi bucket sang trái
            bucketTails[bucketIndex] -= 1

            #showSuffixArray(suffixOffsets)

        # Hậu tố $
        suffixOffsets[0] = len(string)

        #showSuffixArray(suffixOffsets)

        return suffixOffsets

    def makeSuffixArray(self, string, alphabetSize):
        bytestring = string.encode('ascii')
        return self.makeSuffixArrayByInducedSorting(bytestring, alphabetSize)

class BurrowsWheelerTransform:
    k = 100
    c = 100
    d = 0
    def __init__(self, text, k=100, c=100, d=0):
        self.k = k
        self.d = d
        self.c = c
        self.bwt = self.constructBWTfromSA(text)
        self.psa = self.ConstructPartialSuffixarray(text)
        self.alphaMap = self.constructAlphaMap(text)
        self.checkpoints = self.constructCheckpoints()
        self.firstOccurence = self.constructFirstOccurrence()
        
    def constructBWTfromSA(self, Text):
        '''
        Tạo BWT từ mảng hậu tố SA
        '''
        sa = SuffixArray(Text).suffix_array
        bwt = ''
        for i in sa:
            if  i == 0:
                bwt += '$'
            else:
                bwt += Text[i-1]
        return bwt
    
    def ConstructPartialSuffixarray(self, Text):
        '''
        Tạo mảng hậu tố một phần
        '''
        sa = SuffixArray(Text).suffix_array
        partial_suffixarray = {}
        for i in range(len(sa)):
            if sa[i] % self.k == 0:
                partial_suffixarray[i]=sa[i]
        return partial_suffixarray
    
    def constructAlphaMap(self, text):
        '''
        Tạo alphaMap mã hóa các ký tự dưới dạng số từ 0 theo thứ tự từ điển
        '''
        alpha_map = dict()
        alpha_map['$'] = 0
        for i, c in enumerate(sorted(set(text))):
            alpha_map[c] = i+1
        return alpha_map
    
    def constructCheckpoints(self):
        """ 
        Tạo ma trận điểm kiểm tra
        """
        row = len(self.alphaMap) * [0]
        
        counts = [list(row)]  # first copy
        i = 1
        for character in self.bwt:
            
            row[self.alphaMap[character]] += 1  # increase index
            if i % self.c == 0:
                counts.append(list(row))  # append copy
            i += 1
        return counts

    def constructFirstOccurrence(self):
        """ 
        Tạo first occurrence
        """
        firstOccurence = []
        counter = 0
        for character in self.alphaMap.keys():
            firstOccurence.append(counter)
            counter += self.bwt.count(character)
        return firstOccurence

    def BetterBWMatching(self, pattern):
        '''
        Hàm trả về khoảng pattern matching tương ứng với index cột đầu
        '''
        if pattern == '':
            return
        
        top = 0
        bottom = len(self.bwt) - 1
        i = len(pattern)
        while top <= bottom:
            i -= 1
            if i >= 0:
                symbol = pattern[i]
                if symbol not in self.alphaMap.keys():
                    return
                # print(symbol)
                # print(top, bottom)
                top = self.firstOccurence[self.alphaMap[symbol]] + self.count(top, symbol)
                bottom = self.firstOccurence[self.alphaMap[symbol]] + self.count(bottom+1, symbol)-1
                # print(top, bottom)
            else:
                return top, bottom
        
    def findMatchesPosition(self, top, bottom):
        result = []
        
        for i in range(top, bottom+1):
            result.append(self.getSuffixArrayEntry(i))
        return result

    def getSuffixArrayEntry(self, first_index):
        
        counter = 0
        while True:
            
            if first_index in self.psa.keys():
                
                position = self.psa[first_index] + counter
                if position > len(self.bwt) - 2:
                    position = position - len(self.bwt)-1
                return position
            else:
                char = self.bwt[first_index]
                
                first_index = self.firstOccurence[self.alphaMap[char]] + self.count(first_index, char)
                
                counter += 1
       
    def count(self, row, symbol):
        # if row % self.k == 0:
        #     return self.checkpoints[row//self.k][self.alphaMap[symbol]]
        
        checkpoint = row // self.c #hàng checkpoint gần nhất < row
        count = self.checkpoints[checkpoint][self.alphaMap[symbol]]
        for i in range(checkpoint * self.c, row):
            current_symbol = self.bwt[i]
            if current_symbol == symbol:
                count += 1
        return count   

    def approxMatching(self, p1, p2):
        error = 0
        for i in range(len(p1)):
            if p1[i] != p2[i]:
                error += 1
                if error > self.d:
                    return False
        return True

    def findOccurrences(self, patterns):
        #bwt, start = firstO, order = sa, count = checkpoint        
        
        occs = []
        for pattern in patterns:
            currOccs = set()
            n = len(pattern)
            k = n // (self.d+1)
            #seeds: sets(substring, position substring start)
            seeds = [(pattern[i*k:(i+1)*k], i*k) for i in range(self.d)] + [(pattern[self.d*k:n], self.d*k)]
            
            for p, offset in seeds:
                top = 0
                bottom = len(self.bwt) - 1
                currIndex = len(p) - 1
                while top <= bottom:
                    if currIndex >= 0:
                        symbol = p[currIndex]
                        if symbol not in self.alphaMap.keys():
                            break
                        currIndex -= 1
                        if self.count(bottom+1,symbol) - self.count(top, symbol) > 0:
                            top = self.firstOccurence[self.alphaMap[symbol]] + self.count(top, symbol)
                            bottom = self.firstOccurence[self.alphaMap[symbol]] + self.count(bottom+1,symbol) - 1
                        else:
                            break
                    else:
                        for i in range(top, bottom + 1):
                            #currOccs.add(order[i]-offset)
                            currOccs.add(self.getSuffixArrayEntry(i)-offset)
                        break
            for occ in currOccs:
                if self.approxMatching(text[occ:occ+n], pattern):
                    occs.append(occ)
        return occs

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
class Input:
    def __init__(self, genomefile, patternsfile):
        self.refGenome = self.readFasta(genomefile)
        self.patterns = self.readFastq(patternsfile)
        
    def readFasta(self, filename):
        fasta_sequences = SeqIO.parse(open(filename),'fasta') 
        genome_and_seq = []
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            genome_and_seq.append(sequence)
        #Chỉ lấy complete Genome        
        genome = genome_and_seq[0]
        return genome
    
    def readFastq(self, filename):

        patterns = []
        with open(filename) as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                patterns.append(seq)
        return patterns
                
            
if __name__ == '__main__':
    import time
    import sys
    
    # genomefile = "sequence.fasta"
    # patternsfile = "SRR24490774.fastq"
    # input = Input(genomefile, patternsfile)
    # text = input.refGenome
    # patterns = input.patterns
    
    text = "panamabananas"
    patterns = ['nan','ana','ama']
    
    k = 100
    c = k
    
    BWTMatching = BurrowsWheelerTransform(text,k,c)
    
    
    res = []
    for k in range(len(patterns)//5):
        if(BWTMatching.BetterBWMatching(patterns[k])):
            top, bottom = BWTMatching.BetterBWMatching(patterns[k])
            res.append(BWTMatching.findMatchesPosition(top, bottom))
            
    # print(res)

'''
Một vài đoạn trong mã nguồn được viết lại dựa trên sự tham khảo từ các nguồn: 
    https://github.com/rayguo233/read-mapping
    https://zork.net/~st/jottings/sais.html 
    https://github.com/elmar-hinz/Python.Challenges
'''
    
    