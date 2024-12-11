class Matrix:
    def __init__(self, rows):
        self.data = rows  # Initialize with a list of lists
        self.flops = 0
        self.cols = len(rows)
        self.rows = len(rows[0])
    def __getitem__(self, index):
        # Return the row at the given index
        return self.data[index]

    def __setitem__(self, index, value):
        # Set the row at the given index
        self.data[index] = value

    def __repr__(self) -> str:
        output = ''
        for row in self.data:
            output+="|"
            for elem in row:
                output+=' '+"%.16f" %elem+' '
            output+='|\n'
        return output

    def swapItems(self,item1,item2):
        temp = self.data[item1[0]][item1[1]]
        self.data[item1[0]][item1[1]] = self.data[item2[0]][item2[1]]
        self.data[item2[0]][item2[1]] = self.data[item1[0]][item1[1]]

    def swapCols(self,index1,index2):
        for row in self.data:
            temp = row[index1]
            row[index1] = row[index2]
            row[index2] = temp

    def swapRows(self,idx1,idx2):
        temp = self.data[idx1]
        self.data[idx1] = self.data[idx2]
        self.data[idx2] = temp

    def mulRow(self,idx,k):
        self.data[idx] = [elem * k for elem in self.data[idx]]
        
        self.flops += len(self.data[idx])

    def divRow(self,idx,k):
        self.mulRow(idx,1/k)
    
    def mulCol(self,idx,k):
        for row in self.data:
            row[idx] *=k
            self.flops+=1

    def divCol(self,idx,k):
        self.mulCol(idx,1/k)

    def AddRow(self,idx1,idx2):
        for i, item in enumerate(self.data[idx2]):
            self.data[idx1][i] += item
            self.flops+=1

    def SubRow(self,idx1,idx2):
        for i, item in enumerate(self.data[idx2]):
            self.data[idx1][i] -= item
            self.flops+=1

    def gauss(self):
        for i in range(len(self.data)-1):
            for j in range(i+1,len(self.data)):
                self.mulRow(j, self.data[i][i]/self.data[j][i])
                self.SubRow(j,i)

    def PPgauss(self):
        for i in range(len(self.data)-1):
            for j in range(i+1,len(self.data)):
                self.mulRow(j, self.data[i][i]/self.data[j][i])
                self.SubRow(j,i)

    def reverse(self):
        for i in range(len(self.data)-1,0):
            for j in range(i+1,len(self.data)):
                self.mulRow(j, self.data[i][i]/self.data[j][i])
                self.SubRow(j,i)
