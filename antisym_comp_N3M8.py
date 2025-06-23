import numpy as np

from qiskit.circuit import QuantumRegister, ClassicalRegister, AncillaRegister
from qiskit.circuit.library import BlueprintCircuit, C4XGate
from qiskit.circuit.library import MCMT
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Statevector
import qiskit.quantum_info as qi

from qiskit.circuit import QuantumCircuit

class ParallelIntegerComparator(BlueprintCircuit):
    def __init__(
        self,
        size: int,
        geq: bool = True,
        name: str = 'Compare',
    ) -> None:
        """Create a new parallel comparator circuit.
        Args:
            size: Number of qubits required to store each one of the two integers.
            geq: If True, evaluate a ``>=`` condition, otherwise ``<``.
            name: Name of the circuit.
        """
        super().__init__(name=name)
        self._size = None
        self._geq = None
 
        self.size = size
        self.geq = geq
 
    @property
    def size(self) -> int:
        return self._size
 
    @size.setter
    def size(self, size: int) -> None:
        if self._size is None or size != self._size:
            self._invalidate()
            self._size = size
 
            if size is not None:
                # set the new qubit registers
                A = QuantumRegister(size, 'a')
                B = QuantumRegister(size, 'b')
                Anc = AncillaRegister(size, 'anc')
 
                self.qregs = [A, B, Anc]
                num_ancillas = Anc.size
 
    @property
    def geq(self) -> bool:
        return self._geq
 
    @geq.setter
    def geq(self, geq: bool) -> None:
        if geq != self._geq:
            self._invalidate()
            self._geq = geq
 
    @staticmethod
    def _compare2():
        X = QuantumRegister(2, 'X')
        Y = QuantumRegister(2, 'Y')
        W = AncillaRegister(1, 'W')
 
        cmp2 = QuantumCircuit(X, Y, W)
 
        cmp2.x(W)
 
        cmp2.cx(Y[0], X[0])
        cmp2.cx(Y[1], X[1])
        cmp2.cswap(X[1], X[0], W[0])
        cmp2.cswap(X[1], Y[0], Y[1])
        cmp2.cx(Y[0], X[0])
        return cmp2.to_gate(label = 'cmp2')
 
    def _check_configuration(self, raise_on_failure: bool = True) -> bool:
        """Check if the current configuration is valid."""
        valid = True
 
        if self._size is None:
            valid = False
            if raise_on_failure:
                raise AttributeError("size is not set.")
 
        required_num_qubits = 3 * self.size
        if self.num_qubits != required_num_qubits:
            valid = False
            if raise_on_failure:
                raise CircuitError("Number of qubits does not match required number of qubits.")
 
        return valid
 
    def _build(self) -> None:
        """If not already built, build the circuit."""
        if  self._is_built:
            return
 
        super(). _build()
 
        A = self.qregs[0]
        B = self.qregs[1]
        Anc = self.qregs[2]
 
        l = list(range(self.size))
        _sequence = []
        steps = int(np.ceil(np.log2(self.size)))
        anc_idx = 1
        for step in range(steps):
            _sequence.append(l.copy())
            for m in range(int(len(l)/2) -1, -1, -1):
                self.append(self._compare2(), [A[l[2 * m]], A[l[2 * m + 1]], B[l[2 * m]], B[l[2 * m + 1]], Anc[anc_idx]])
                l.pop(2 * m + 1)
                anc_idx = anc_idx + 1
 
        # Copy result
        self.ccx(A[0], B[0], Anc[0], ctrl_state='10')
        if self._geq:
            self.x(Anc[0])
 
        # Uncompute
        for step in range(steps):
            l = _sequence.pop(-1)
            for m in range(int(len(l)/2)):
                anc_idx = anc_idx - 1
                self.append(self._compare2().inverse(), [A[l[2 * m]], A[l[2 * m + 1]], B[l[2 * m]], B[l[2 * m + 1]], Anc[anc_idx]])

def controlled_swap(qc,n1,n2,nq,jj_ancilla,n0_ancilla):
    # swap particles
    n_ancilla=len(jj_ancilla)
    for i in range(n_ancilla):
        if jj_ancilla[i]==0:
            qc.x(n0_ancilla+i)
    anc=list(range(n0_ancilla,n0_ancilla+n_ancilla))
    if n_ancilla>1:
        pos_anc=n0_ancilla+n_ancilla
        qc.mcx(anc,pos_anc)
    else:
        pos_anc=n0_ancilla
    for i in range(nq):
        qc.cswap(pos_anc,i+n1*nq,i+n2*nq)
    if n_ancilla>1:
        qc.mcx(anc,pos_anc)
    for i in range(n_ancilla):
        if jj_ancilla[i]==0:
            qc.x(n0_ancilla+i)
    return qc


def uniform_superposition_qubit_states(qc,nq, qlist, M):
    # from reference arXiv:2306.11747v2
    if M>2**nq:
        print("Please adjust the number of qubits available")
        return

    if nq > len(qlist):
        print("Please adjust the number of qubits available")
        return
        
    m_=bin(M)[2:]
    l_list=[]
    ipos=0
    #print(m_)
    for c in m_[::-1]:
        if c=="1":
            #print("ipos: ", ipos)
            #print("qlist[ipos]: ", qlist[ipos])
            l_list = l_list + [ipos]
        ipos+=1
    for l in l_list[1:]:
        qc.x(qlist[l])

    M0=2**l_list[0]
    for l in range(l_list[0]):
        qc.h(qlist[l])

    theta=-2.*np.arccos(np.sqrt(M0/M))
    #print("l_list = ", l_list)
    if len(l_list) >= 2:
        qc.ry(theta,qlist[l_list[1]])
        qc.x(qlist[l_list[1]])
        for l in range(l_list[0],l_list[1]):
            qc.ch(qlist[l_list[1]],qlist[l])
        qc.x(qlist[l_list[1]])
    #print("loop-range - ", len(l_list)-1)
    #print("List of predicted iters = ", [i for i in range(1,len(l_list)-1)])
    for m in range(1,len(l_list)-1):
        theta=-2.*np.arccos(np.sqrt((2**l_list[m])/(M-M0)))
        qc.x(qlist[l_list[m]])
        qc.cry(theta,qlist[l_list[m]],qlist[l_list[m+1]])
        qc.x(qlist[l_list[m]])
        qc.x(qlist[l_list[m+1]])
        for l in range(l_list[m],l_list[m+1]):
            qc.ch(qlist[l_list[m+1]],qlist[l])
        qc.x(qlist[l_list[m+1]])
        M0+=2**l_list[m]
                       
    return qc

# The conjugate-transpose of the above
def uniform_superposition_qubit_states_uncompute(qc,nq, qlist, M):
    if M>2**nq:
        print("Please adjust the number of qubits available")
        return

    if nq > len(qlist):
        print("Please adjust the number of qubits available")
        return

    m_=bin(M)[2:]
    l_list=[]
    ipos=0
    #print(m_)
    for c in m_[::-1]:
        if c=="1":
            #print("ipos: ", ipos)
            #print("qlist[ipos]: ", qlist[ipos])
            l_list = l_list + [ipos]
        ipos+=1

    M0=2**l_list[0]

    #print("M0 = ", M0)

    #print("loop-range - ", len(l_list)-1)
    
    for m in range(1,len(l_list)-1):
        #print("m = ", m)
        M0+=2**l_list[m]

    #print("M0 = ", M0)
    
    for m in range(len(l_list)-2,0,-1):
        #print("m = ", m)
        M0-=2**l_list[m]
        theta=-2.*np.arccos(np.sqrt((2**l_list[m])/(M-M0)))
        qc.x(qlist[l_list[m+1]])
        for l in range(l_list[m+1]-1,l_list[m]-1,-1):
            qc.ch(qlist[l_list[m+1]],qlist[l])
        qc.x(qlist[l_list[m+1]])
        qc.x(qlist[l_list[m]])
        qc.cry(-theta,qlist[l_list[m]],qlist[l_list[m+1]])
        qc.x(qlist[l_list[m]])

    theta=-2.*np.arccos(np.sqrt(M0/M))
    if len(l_list) >= 2:
        qc.x(qlist[l_list[1]])
        for l in range(l_list[1]-1,l_list[0]-1,-1):
            qc.ch(qlist[l_list[1]],qlist[l])
        qc.x(qlist[l_list[1]])
        qc.ry(-theta,qlist[l_list[1]])
            
    for l in range(l_list[0]):
        qc.h(qlist[l])

    for l in l_list[1:]:
        qc.x(qlist[l])

    return qc

def add_dict(d1,d2):
    d={}
    for k in d1.keys():
        d[k]=d1[k]
    for k in d2.keys():
        if k in d.keys():
            d[k]+=d2[k]
        else:
            d[k]=d2[k]
    return d

def sub_dict(d1,d2):
    for k in d2.keys():
        d2[k]=-d2[k]
    d=add_dict(d1,d2)
    for k in d2.keys():
        d2[k]=-d2[k]
    return d
    

def mul_dict(d1,d2,remove_collisions=False):
    d={}
    for k1 in d1.keys():
        for k2 in d2.keys():
            if remove_collisions and k1==k2:
                continue
            d[k1+k2]=d1[k1]*d2[k2]
    return d
    
def circ_phi1(qc,theta,n_shift=0):
    qc.ry(theta,n_shift+1)
    qc.x(n_shift+1)
    qc.cx(n_shift+1,n_shift)
    qc.x(n_shift+1)
    return qc

def ccirc_phi1(qc,theta,nc,n_shift=0):
    qc.cry(theta,nc,n_shift+1)
    qc.cx(nc,n_shift+1)
    qc.ccx(nc,n_shift+1,n_shift)
    qc.cx(nc,n_shift+1)
    return qc

def ccirc_phi1_inv(qc,theta,nc,n_shift=0):
    qc.cx(nc,n_shift+1)
    qc.ccx(nc,n_shift+1,n_shift)
    qc.cx(nc,n_shift+1)
    qc.cry(-theta,nc,n_shift+1)
    return qc


def circ_phi1_inv(qc,theta,n_shift=0):
    qc.x(n_shift+1)
    qc.cx(n_shift+1,n_shift)
    qc.x(n_shift+1)
    qc.ry(-theta,n_shift+1)
    return qc

def circ_phi2(qc,theta_0,theta_1,n0,n_shift=0):
    qc.x(n0+n_shift)
    qc.ry(theta_0,n_shift)
    qc.x(n_shift)
    qc.cry(theta_1,n_shift,n_shift+1)
    qc.x(n_shift)
    return qc

def circ_phi2_inv(qc,theta_0,theta_1,n0,n_shift=0):
    qc.x(n0+n_shift)
    qc.x(n_shift)
    qc.cry(-theta_1,n_shift,n_shift+1)
    qc.x(n_shift)
    qc.ry(-theta_0,n_shift)
    return qc

# -------------------------------------------
# Evan's code for Abrams & Lloyd algorithm
# -------------------------------------------
n = 3 # particles
nm = 3 # qubits per particle
size = n*nm

# Number of qubits in register A
nQA = size
# Number of qubits in register B
nQB = 6
# Number of qubits in register C
nQC = 6

# Number of ancilla qubits to store the result of sorting
nQSA = 3
# Number of "scratch" qubits required for integer comparison but don't store any results
nQSS = 2

# Total number of qubits
nQ = nQA + nQB + nQC + nQSA + nQSS

# Compare 2 integers stored on 3 qubits each
comparatorA = ParallelIntegerComparator(nm)

# B and C registers only need 2 qubits per particle to encode permutations of 1,2,3
comparatorBC = ParallelIntegerComparator(2)

circ = QuantumCircuit(nQ, nQ - 1)
print("Num total Qubits: ", nQ)

# Prepare A in |psi> = |0 1 2>
circ.x(3)
circ.x(7)

# Prepare B in |0 1 2> + |0 2 1> + |1 0 2> + |1 2 0> + |2 0 1> + |2 1 0>
circ = uniform_superposition_qubit_states(circ, 2, [9,10], 3)
circ = uniform_superposition_qubit_states(circ, 2, [11,12], 2)

c4x1 = MCMT('x', num_ctrl_qubits=4, num_target_qubits=1)
c4x2 = MCMT('x', num_ctrl_qubits=4, num_target_qubits=2)

circ.x(9)
circ.x(10)
circ.x(11)
circ.x(12)
circ.append(c4x1, [9,10,11,12,14])
circ.x(9)
circ.x(10)
circ.x(11)
circ.x(12)

circ.x(9)
circ.x(10)
circ.x(13)
circ.append(c4x1, [9,10,13,14,11])
circ.x(9)
circ.x(10)
circ.x(13)

circ.x(9)
circ.x(10)
circ.x(13)
circ.x(14)
circ.append(c4x2, [9,10,13,14,11,12])
circ.x(9)
circ.x(10)
circ.x(13)
circ.x(14)

circ.x(9)
circ.x(10)
circ.x(11)
circ.append(c4x1, [9,10,11,12,13])
circ.x(9)
circ.x(10)
circ.x(11)

circ.x(9)
circ.x(11)
circ.x(12)
circ.append(c4x1, [9,10,11,12,13])
circ.x(9)
circ.x(11)
circ.x(12)

circ.x(10)
circ.x(11)
circ.x(12)
circ.append(c4x1, [9,10,11,12,14])
circ.x(10)
circ.x(11)
circ.x(12)

circ.x(10)
circ.x(13)
circ.x(14)
circ.append(c4x2, [9,10,13,14,11,12])
circ.x(10)
circ.x(13)
circ.x(14)

# Prepare C in |0 1 2>
circ.x(17)
circ.x(20)

# Sort B
circ.append(comparatorBC, [9,10,11,12,21,24])
circ.cswap(21, 9, 11)
circ.cswap(21, 10, 12)

circ.append(comparatorBC, [11,12,13,14,22,24])
circ.cswap(22, 11, 13)
circ.cswap(22, 12, 14)

circ.append(comparatorBC, [9,10,11,12,23,24])
circ.cswap(23, 9, 11)
circ.cswap(23, 10, 12)

# Apply phase to swapped components
circ.z(21)
circ.z(22)
circ.z(23)

# Apply same transformation to A and C
circ.cswap(21, 0, 3)
circ.cswap(21, 1, 4)
circ.cswap(21, 2, 5)
circ.cswap(22, 3, 6)
circ.cswap(22, 4, 7)
circ.cswap(22, 5, 8)
circ.cswap(23, 0, 3)
circ.cswap(23, 1, 4)
circ.cswap(23, 2, 5)

circ.cswap(21, 15, 17)
circ.cswap(21, 16, 18)
circ.cswap(22, 17, 19)
circ.cswap(22, 18, 20)
circ.cswap(23, 15, 17)
circ.cswap(23, 16, 18)

# Reverse sort B
circ.cswap(23, 9, 11)
circ.cswap(23, 10, 12)
circ.append(comparatorBC, [9,10,11,12,23,24])

circ.cswap(22, 11, 13)
circ.cswap(22, 12, 14)
circ.append(comparatorBC, [11,12,13,14,22,24])

circ.cswap(21, 9, 11)
circ.cswap(21, 10, 12)
circ.append(comparatorBC, [9,10,11,12,21,24])

# Zero out B by controlling on C
circ.x(17)
circ.x(18)
circ.ccx(17, 18, 9)
circ.x(18)
circ.ccx(17,18,13)
circ.x(17)
circ.x(18)
circ.ccx(17, 18, 11)
circ.x(18)

circ.x(19)
circ.x(20)
circ.ccx(19, 20, 10)
circ.x(20)
circ.ccx(19,20,14)
circ.x(19)
circ.x(20)
circ.ccx(19, 20, 12)
circ.x(20)

# Sort A and B
circ.append(comparatorBC, [15,16,17,18,21,24])
circ.cswap(21, 15, 17)
circ.cswap(21, 16, 18)

circ.cswap(21, 0, 3)
circ.cswap(21, 1, 4)
circ.cswap(21, 2, 5)

circ.append(comparatorBC, [17,18,19,20,22,24])
circ.cswap(22, 17, 19)
circ.cswap(22, 18, 20)

circ.cswap(22, 3, 6)
circ.cswap(22, 4, 7)
circ.cswap(22, 5, 8)

circ.append(comparatorBC, [15,16,17,18,23,24])
circ.cswap(23, 15, 17)
circ.cswap(23, 16, 18)

circ.cswap(23, 0, 3)
circ.cswap(23, 1, 4)
circ.cswap(23, 2, 5)

# Zero out C
circ.x(17)
circ.x(20)

# Reverse sort A
circ.cswap(23, 0, 3)
circ.cswap(23, 1, 4)
circ.cswap(23, 2, 5)
circ.append(comparatorA, [0,1,2,3,4,5,23,24,25])

circ.cswap(22, 3, 6)
circ.cswap(22, 4, 7)
circ.cswap(22, 5, 8)
circ.append(comparatorA, [3,4,5,6,7,8,22,24,25])

circ.cswap(21, 0, 3)
circ.cswap(21, 1, 4)
circ.cswap(21, 2, 5)
circ.append(comparatorA, [0,1,2,3,4,5,21,24,25])

# Hadamard test for phase under particle exchanges
circ.h(24)
#circ.cswap(24, 0, 3)
#circ.cswap(24, 1, 4)
#circ.cswap(24, 2, 5)

circ.cswap(24, 3, 6)
circ.cswap(24, 4, 7)
circ.cswap(24, 5, 8)


circ.h(24)

circ.measure(range(nQ - 1), range(nQ - 1))
simulator = AerSimulator()
circ = transpile(circ, simulator, basis_gates=["cx",'x','z',"h","rz"], optimization_level=2)
print(dict(circ.count_ops()))
result = simulator.run(circ).result()
counts = result.get_counts(circ)
print(counts)
print("A1  A2  A3  B1  B2  B3  C1  C2  C3  SA1 SA2 SA3 SS")
for bitstring in sorted(counts):
    formatted_bitstring = '   '.join([
        str(int(bitstring[22:25], 2)),
        str(int(bitstring[19:22], 2)),
        str(int(bitstring[16:19], 2)),
        
        str(int(bitstring[14:16], 2)),
        str(int(bitstring[12:14], 2)),
        str(int(bitstring[10:12], 2)),
        
        str(int(bitstring[8:10], 2)),
        str(int(bitstring[6:8], 2)),
        str(int(bitstring[4:6], 2)),
        
        bitstring[3],
        bitstring[2],
        bitstring[1],
        bitstring[0]
    ])
    print(formatted_bitstring)

#-----------------------------------
# Ionel's code for his algorithm
#-----------------------------------
Nq1=3
A=3
sp_states=[[0,0,0],[1,0,0],[0,1,0]]

qreg=QuantumRegister(A*Nq1,"q")
qanc=QuantumRegister(3,'a')
qclas=ClassicalRegister(10,'c')

qc=QuantumCircuit(qreg,qanc, qclas)
qc.h(A*Nq1)
qc.z(A*Nq1)
for n in range(A):
    for i in range(Nq1):
        if sp_states[n][i]==1:
            qc.x(n*Nq1+i)

#for i in range(Nq1):
#    qc.cswap(A*Nq1,i,i+Nq1)

qc=controlled_swap(qc,0,1,Nq1,[1],A*Nq1)


for i in range(Nq1):
    if sp_states[1][i]==0:
        qc.x(i)
    
qc.mcx(list(range(Nq1)),A*Nq1)
for i in range(Nq1):
    if sp_states[1][i]==0:
        qc.x(i)

#qc.barrier()
#three particles now

qc=uniform_superposition_qubit_states(qc,2, [A*Nq1,A*Nq1+1], 3)
qc.z(A*Nq1)
qc.z(A*Nq1+1)


qc=controlled_swap(qc,0,2,Nq1,[0,1],A*Nq1)
qc=controlled_swap(qc,1,2,Nq1,[1,0],A*Nq1)


for i in range(Nq1):
#    qc.x(i)
#    qc.x(i+Nq1)
    if sp_states[2][i]==0:
        qc.x(i)
        qc.x(i+Nq1)

qc.mcx(list(range(Nq1,2*Nq1)),A*Nq1)
qc.mcx(list(range(Nq1)),A*Nq1+1)

for i in range(Nq1):
#    qc.x(i)
#    qc.x(i+Nq1)
    if sp_states[2][i]==0:
        qc.x(i)
        qc.x(i+Nq1)

print("\n")

# Hadamard test for phase under particle exchanges
qc.h(9)

qc.cswap(9, 3, 6)
qc.cswap(9, 4, 7)
qc.cswap(9, 5, 8)

qc.h(9)

qc.measure(range(10), range(10))
simulator = AerSimulator()
qc = transpile(qc, simulator, basis_gates=["cx",'x','z',"h","rz"], optimization_level=2)
print(dict(qc.count_ops()))
result = simulator.run(qc).result()
counts = result.get_counts(qc)
print(counts)
print("P1  P2  P3  ANC")
for bitstring in sorted(counts):
    formatted_bitstring = '   '.join([
        str(int(bitstring[7:10], 2)),
        str(int(bitstring[4:7], 2)),
        str(int(bitstring[1:4], 2)),
        str(int(bitstring[0], 2))
    ])
    print(formatted_bitstring)

