import numpy as np

from qiskit.circuit import QuantumRegister, ClassicalRegister, AncillaRegister
from qiskit.circuit.library import BlueprintCircuit, C4XGate
from qiskit.circuit.library import MCMT
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Statevector

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

n = 2 # particles
nm = 2 # qubits per particle, 4 single-particle states
size = n*nm

comparator1 = ParallelIntegerComparator(1)
comparator = ParallelIntegerComparator(nm)
num_comp_anc = comparator.num_qubits - size
print("Num Extra Comparator Qubits: ", num_comp_anc) 

q_a1 = QuantumRegister(2, name="P1")
q_a2 = QuantumRegister(2, name="P2")
q_b = QuantumRegister(2, name="B")
q_c = QuantumRegister(2, name="C")
q_anc = QuantumRegister(2,name="Anc")
c_a = ClassicalRegister(10)

circ = QuantumCircuit(q_a1, q_a2, q_b, q_c, q_anc, c_a)
#circ = QuantumCircuit(3*(size + 1) + num_comp_anc, 3*(size + 1))
print("Num total Qubits: ", 3*(size + 1) + num_comp_anc - 1)
#circ.h(range(2 * size))

# Prepare A in |psi> = |1 2>
#circ.x(0)
circ.x(2)
#circ.x(3)
#circ.x(10)
#circ.x(11)

# Prepare B in 1/2*( |1 2> + |2 1> )
circ.h(4)
circ.x(4)
circ.cx(4, 5)
circ.x(4)

# Prepare C in |1 2>
circ.x(7)

# Sort B
circ.append(comparator1, [4,5,8])
circ.cswap(8, 4, 5)

# Apply phase to swapped component
circ.z(8)

# Perform same operations on A and C
circ.cswap(8, 0, 2)
circ.cswap(8, 1, 3)
circ.cswap(8, 6, 7)

# Reverse sort on B
circ.cswap(8, 4, 5)

# Compare to reset B ancilla
circ.append(comparator1, [4,5,8])

# Zero out B by controlling on C
circ.cx(6, 4)
circ.cx(7, 5)

# Sort both A and C
circ.append(comparator1, [6,7,8])
circ.cswap(8, 0, 2)
circ.cswap(8, 1, 3)

circ.cswap(8, 6, 7)

# Zero out C registers using fact that C is in |1 2>
circ.x(7)

# Reverse sort on A
circ.cswap(8, 0, 2)
circ.cswap(8, 1, 3)
circ.append(comparator, [0,1,2,3,8,9])

print(circ)

# Hadamard test for antisymmetry
circ.h(9)
circ.cswap(9, 0, 2)
circ.cswap(9, 1, 3)
circ.h(9)

#statevector = Statevector(circ)
#state_dict = statevector.to_dict()
#print(state_dict)
#for i in range(nm):
#    circ.cswap(size, i, i + nm)
circ.measure(range(10), range(10))
simulator = AerSimulator()
circ = transpile(circ, simulator)
print(dict(circ.count_ops()))
result = simulator.run(circ).result()
counts = result.get_counts(circ)
print(counts)
print("A1  A2  B1  B2  C1  C2  SA1 SS")
for bitstring in sorted(counts):
    formatted_bitstring = '   '.join([
        str(int(bitstring[8:10], 2)),
        str(int(bitstring[6:8], 2)),
        bitstring[5],
        bitstring[4],
        bitstring[3],
        bitstring[2],
        bitstring[1],
        bitstring[0]
    ])
    print(formatted_bitstring)
#circ.draw('mpl', filename="compcirc.png")

print("Ionel's algorithm: ")
print("-------------------------------------")
circ = QuantumCircuit(5, 5)
circ.h(4)
circ.z(4)
circ.x(2)
circ.cswap(4, 0, 2)
circ.cswap(4, 1, 3)
circ.x(0)

circ.x(0)
circ.cx(0, 4)
circ.x(0)

circ.x(0)

circ.h(4)
circ.cswap(4, 0, 2)
circ.cswap(4, 1, 3)
circ.h(4)


circ.measure(range(5), range(5))
simulator = AerSimulator()
circ = transpile(circ, simulator)
print(dict(circ.count_ops()))
result = simulator.run(circ).result()
counts = result.get_counts(circ)
print(counts)
print("A1  A2 SS")
for bitstring in sorted(counts):
    formatted_bitstring = '   '.join([
        str(int(bitstring[3:5], 2)),
        str(int(bitstring[1:3], 2)),
        bitstring[0]
    ])
    print(formatted_bitstring)


n = 3 # particles
nm = 2 # qubits per particle
size = n*nm

comparator = ParallelIntegerComparator(nm)
num_comp_anc = comparator.num_qubits - 2*nm
print("Num Extra Comparator Qubits: ", num_comp_anc) 

circ = QuantumCircuit(3*(size) + 3 + num_comp_anc, 3*(size) + 3 + 1)
print("Num total Qubits: ", 3*(size) + 3 + num_comp_anc - 1)
#circ.h(range(2 * size))

# Prepare A in |psi> = |0 1 2>
circ.x(2)
circ.x(5)

# Prepare B in |1 2 3> + |1 3 2> + |2 1 3> + |2 3 1> + |3 1 2> + |3 2 1>

circ.h(8)

m = 3
l0 = 0
l1 = 1
m0 = 2**l0
circ.x(6)
circ.ry(-2*np.arccos(np.sqrt(m0/m)), 6)
circ.x(6)
circ.ch(6, 7)
circ.x(6)


c4x1 = MCMT('x', num_ctrl_qubits=4, num_target_qubits=1)
c4x2 = MCMT('x', num_ctrl_qubits=4, num_target_qubits=2)

circ.x(6)
circ.x(7)
circ.x(8)
circ.x(9)
circ.append(c4x1, [6,7,8,9,11])
circ.x(6)
circ.x(7)
circ.x(8)
circ.x(9)

circ.x(6)
circ.x(7)
circ.x(10)
circ.append(c4x1, [6,7,10,11,8])
circ.x(6)
circ.x(7)
circ.x(10)

circ.x(6)
circ.x(7)
circ.x(10)
circ.x(11)
circ.append(c4x2, [6,7,10,11,8,9])
circ.x(6)
circ.x(7)
circ.x(10)
circ.x(11)

circ.x(6)
circ.x(7)
circ.x(8)
circ.append(c4x1, [6,7,8,9,10])
circ.x(6)
circ.x(7)
circ.x(8)

circ.x(6)
circ.x(8)
circ.x(9)
circ.append(c4x1, [6,7,8,9,10])
circ.x(6)
circ.x(8)
circ.x(9)

circ.x(7)
circ.x(8)
circ.x(9)
circ.append(c4x1, [6,7,8,9,11])
circ.x(7)
circ.x(8)
circ.x(9)

circ.x(7)
circ.x(10)
circ.x(11)
circ.append(c4x2, [6,7,10,11,8,9])
circ.x(7)
circ.x(10)
circ.x(11)

# Prepare C in |1 2 3>
circ.x(14)
circ.x(17)

# Sort B
circ.append(comparator, [6,7,8,9,18,21])
circ.cswap(18, 6, 8)
circ.cswap(18, 7, 9)

circ.append(comparator, [8,9,10,11,19,21])
circ.cswap(19, 8, 10)
circ.cswap(19, 9, 11)

circ.append(comparator, [6,7,8,9,20,21])
circ.cswap(20, 6, 8)
circ.cswap(20, 7, 9)

# Apply phase to swapped components
circ.z(18)
circ.z(19)
circ.z(20)

# Apply same transformation to A and C
circ.cswap(18, 0, 2)
circ.cswap(18, 1, 3)
circ.cswap(19, 2, 4)
circ.cswap(19, 3, 5)
circ.cswap(20, 0, 2)
circ.cswap(20, 1, 3)

circ.cswap(18, 12, 14)
circ.cswap(18, 13, 15)
circ.cswap(19, 14, 16)
circ.cswap(19, 15, 17)
circ.cswap(20, 12, 14)
circ.cswap(20, 13, 15)

# Reverse sort B
circ.cswap(20, 6, 8)
circ.cswap(20, 7, 9)
circ.append(comparator, [6,7,8,9,20,21])

circ.cswap(19, 8, 10)
circ.cswap(19, 9, 11)
circ.append(comparator, [8,9,10,11,19,21])

circ.cswap(18, 6, 8)
circ.cswap(18, 7, 9)
circ.append(comparator, [6,7,8,9,18,21])

# Zero out B by controlling on C
circ.x(14)
circ.x(15)
circ.ccx(14, 15, 6)
circ.x(15)
circ.ccx(14,15,10)
circ.x(14)
circ.x(15)
circ.ccx(14, 15, 8)
circ.x(15)

circ.x(16)
circ.x(17)
circ.ccx(16, 17, 7)
circ.x(17)
circ.ccx(16,17,11)
circ.x(16)
circ.x(17)
circ.ccx(16, 17, 9)
circ.x(17)

# Sort A and B
circ.append(comparator, [0,1,2,3,18,21])
circ.cswap(18, 0, 2)
circ.cswap(18, 1, 3)

circ.cswap(18, 12, 14)
circ.cswap(18, 13, 15)

circ.append(comparator, [2,3,4,5,19,21])
circ.cswap(19, 2, 4)
circ.cswap(19, 3, 5)

circ.cswap(19, 14, 16)
circ.cswap(19, 15, 17)

circ.append(comparator, [0,1,2,3,20,21])
circ.cswap(20, 0, 2)
circ.cswap(20, 1, 3)

circ.cswap(20, 12, 14)
circ.cswap(20, 13, 15)

# Zero out C
circ.x(14)
circ.x(17)

# Reverse sort A
circ.cswap(20, 0, 2)
circ.cswap(20, 1, 3)
circ.append(comparator, [0,1,2,3,20,21])

circ.cswap(19, 2, 4)
circ.cswap(19, 3, 5)
circ.append(comparator, [2,3,4,5,19,21])

circ.cswap(18, 0, 2)
circ.cswap(18, 1, 3)
circ.append(comparator, [0,1,2,3,18,21])

# Hadamard test for phase under particle exchanges
circ.h(21)
circ.cswap(21, 2, 4)
circ.cswap(21, 3, 5)

circ.h(21)

#statevector = Statevector(circ)
#state_dict = statevector.to_dict(decimals=4)
#print(state_dict)
#for i in range(nm):
#    circ.cswap(size, i, i + nm)
circ.measure(range(3*(size) + 3 + 1), range(3*(size) + 3 + 1))
simulator = AerSimulator()
circ = transpile(circ, simulator)
#print(circ)
print(dict(circ.count_ops()))
result = simulator.run(circ).result()
counts = result.get_counts(circ)
print(counts)
print("A1  A2  A3  B1  B2  B3  C1  C2  C3  SA1 SA2 SA3 SS")
for bitstring in sorted(counts):
    formatted_bitstring = '   '.join([
        str(int(bitstring[20:22], 2)),
        str(int(bitstring[18:20], 2)),
        str(int(bitstring[16:18], 2)),
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


n = 3 # particles
nm = 3 # qubits per particle
size = n*nm

nQA = size
nQB = 6
nQC = 6

nQSA = 3
nQSS = 2

nQ = nQA + nQB + nQC + nQSA + nQSS

comparatorA = ParallelIntegerComparator(nm)
comparatorBC = ParallelIntegerComparator(2)

num_comp_anc = comparator.num_qubits - 2*nm
print("Num Extra Comparator Qubits: ", num_comp_anc) 

circ = QuantumCircuit(nQ, nQ - 1)
print("Num total Qubits: ", nQ)
#circ.h(range(2 * size))

# Prepare A in |psi> = |0 1 2>
circ.x(3)
circ.x(7)

# Prepare B in |1 2 3> + |1 3 2> + |2 1 3> + |2 3 1> + |3 1 2> + |3 2 1>
circ = uniform_superposition_qubit_states(circ, 2, [9,10], 3)
circ = uniform_superposition_qubit_states(circ, 2, [11,12], 2)

#circ.h(11)

m = 3
l0 = 0
l1 = 1
m0 = 2**l0
#circ.x(9)
#circ.ry(-2*np.arccos(np.sqrt(m0/m)), 9)
#circ.x(10)
#circ.ch(9, 10)
#circ.x(9)


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

# Prepare C in |1 2 3>
circ.x(17)
circ.x(20)

# Sort B
circ.append(comparatorBC, [9,10,11,12,21,24])
circ.cswap(21, 9, 11)
circ.cswap(21, 10, 12)

circ.append(comparatorBC, [11,12,13,14,22,24])
circ.cswap(22, 11, 13)
circ.cswap(22, 12, 14)

circ.append(comparator, [9,10,11,12,23,24])
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

#statevector = Statevector(circ)
#state_dict = statevector.to_dict(decimals=4)
#print(state_dict)
#for i in range(nm):
#    circ.cswap(size, i, i + nm)
circ.measure(range(nQ - 1), range(nQ - 1))
simulator = AerSimulator()
circ = transpile(circ, simulator)
#print(circ)
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

