.data
A:      .float 1,1,1,1,1,10         # Example matrix A (row-major format)
        .float 0,5,3,6,2,13
        .float 0,4,6,7,5,3
        .float 0,0,0,0,5,7
        .float 0,0,0,0,9,8


n:      .word 5                            # Number of equations
m:      .word 6                                 #Augmented matrix (#cols)

infinite_string: .string "Infinitely many solution exist\n"
uniquesol_string: .string "Unique Solutions, check values memory (0x10000080 to 0x10000090)\n"
nosol_string: .string "No solution exists\n"

.text
.globl main
main:
    # Load n into a register
    lw a0, n

    #Load m
    lw s0, m

    # Load the address of matrix A into a2
    la a2, A

    # Load the address of matrix B into a3
    #la a3, B
    
	# Load the address of pivot element  A[i][i] in a4
    mv a4, a2

	#Loop over rows of A matrix
	li a1, 0
    
    # 4 row computations
    li t3, 4

    # Loop over rows
    loop_rows:
        #If a1>=COLS-1
        bge a1, t3, loop_end_echelon_complete

        # Loop over columns for pivot
        li a5, 0

        # Load A[i][i]
        flw fa0, 0(a4)
        
        loop_cols:
            # Load A[i][j]
            flw fa1, 0(a2)

            #This means that is an element is 0, then continue to next column
            feq.s t1, f0, fa0
            bnez t1, loop_cols_end

            #Check for fa1 also, if it is zero
            feq.s t1, f0, fa1
            bnez t1, loop_cols_end

			# Calculate the ratio for pivot (A[i][j] / A[i][i])
            fdiv.s fa2, fa1, fa0

            # Store it back in the memory location
            fsw fa2, 0(a2)

			loop_cols_end:
                # Move to the next column
                addi a2, a2, 4
                addi a5, a5, 1
                bgt s0, a5, loop_cols

        #Here we have to eliminate other rows

        #Memory location for k->a2 and for i we define
        mv t2, a2

        #Loop counter for k
        mv a6, a1

        loop_k:
            bge a6, a0, next_row #next row computation
            beq a1, a6, next_k  #if k and i are equal, continue
            
            #Loop j for traversing and eliminating
            li a7, 0
            
            loop_back_from_zerofactor:
                #Check if 0(t2) is zero or not, if zero next element in column
                
                #Calculate factor (equal to A[k][i] since i have already normalized A[i] row)
                flw fa6, 0(t2)
                feq.s t1, fa6, f1
                bnez t1, zero_factor

                loop_j:
                    bge a7, s0, next_k

                    #Load A[k][j]
                    flw fa7, 0(t2)

                    #Load A[i][j]
                    flw f18, 0(a2)

                    #now temp = factor*A[i][j]
                    fmul.s f19, fa6, f18

                    # A[k][j] - temp
                    fsub.s f19, fa7, f19

                    #Store back to Memory
                    fsw f19, 0(t2)

                    #increment Memory
                    addi t2, t2, 4
                    addi a2, a2, 4
                    addi a7, a7, 1
                    j loop_j

                next_k:
                    addi a6, a6, 1
                    addi a2, a2, -24
                    j loop_k
                
                zero_factor:
                    addi t2, t2, 4
                    addi a2, a2, 4
                    addi a7, a7, 1
                    j loop_back_from_zerofactor

                next_row:

        # Move to the next row
        addi a4, a4, 28
		addi a1, a1, 1

        #De-normalize the ith row (earlier was divided by A[i][i])
        li t2, 0
        loop_denormalize:
            bge t2, s0, denormal_complete

            # fa0 conatins the A[i][i], so multiply it with each element 
            flw fa7, 0(a2)
            
            #If 0(a2) has 0, then skip the multiplication
            feq.s t1, f0, fa7
            bnez t1, loop_skip

            #Multiply
            fmul.s fa7, fa7, fa0
            
            # Store back in memory
            fsw fa7, 0(a2)

            loop_skip:
                #Increment counters and memory
                addi t2, t2, 1
                addi a2, a2, 4
                j loop_denormalize

        denormal_complete:
            bgt a0, a1, loop_rows 
	
    loop_end_echelon_complete:
        #Load address of matrix A
        la a2, A

        #Lets check last row of echelon matrix

        #load A[5][5]
        flw fa0, 112(a2)
        flw f9, 116(a2)
        feq.s t1, fa0, f0
        bnez t1, c1              # branches if A[5][5] =0

        check1: 
            #1st solution
            fdiv.s fa0, f9, fa0
            fsw fa0, 128(a2)

        #Solve for next
        flw f9, 88(a2)
        fmul.s f9, fa0, f9
        flw fa1, 92(a2)
        fsub.s fa1, fa1, f9
        flw f9, 84(a2)
        feq.s t1, f9, f0
        bnez t1, c2

        check2:
            #2nd solution
            fdiv.s fa1, fa1, f9
            fsw fa1, 132(a2)

        #Solve for next
        flw f9, 64(a2)
        fmul.s f9, fa0, f9
        flw fa2, 68(a2)
        fsub.s fa2, fa2, f9
        flw f9, 60(a2)
        fmul.s f9, f9, fa1
        fsub.s fa2, fa2, f9
        flw f9, 56(a2)
        feq.s t1, f9, f0
        bnez t1, c3

        check3:
            #3rd solution
            fdiv.s fa2, fa2, f9
            fsw fa2, 136(a2)

        #Solve for next
        flw f9, 40(a2)
        fmul.s f9, fa0, f9
        flw fa3, 44(a2)
        fsub.s fa3, fa3, f9
        flw f9, 36(a2)
        fmul.s f9, f9, fa1
        fsub.s fa3, fa3, f9
        flw f9, 32(a2)
        fmul.s f9, f9, fa2
        fsub.s fa3, fa3, f9
        flw f9, 28(a2)
        feq.s t1, f9, f0
        bnez t1, c4

        check4:
            #4th solution
            fdiv.s fa3, fa3, f9
            fsw fa3, 140(a2)

        #Solve for next
        flw f9, 16(a2)
        fmul.s f9, fa0, f9
        flw fa4, 20(a2)
        fsub.s fa4, fa4, f9
        flw f9, 12(a2)
        fmul.s f9, f9, fa1
        fsub.s fa4, fa4, f9
        flw f9, 8(a2)
        fmul.s f9, f9, fa2
        fsub.s fa4, fa4, f9
        flw f9, 4(a2)
        fmul.s f9, f9, fa3
        fsub.s fa4, fa4, f9
        flw f9, 0(a2)
        feq.s t1, f9, f0
        bnez t1, c5

        check5:
            #5th solution
            fdiv.s fa4, fa4, f9
            fsw fa4, 144(a2)

        li a0, 4
        la a1, uniquesol_string
        ecall

        j end_program

        c1:
            feq.s t1, f9, f0
            bnez t1, infinite_sols               # if true, infinite solutions
            beqz t1, no_sols                    #if true, no solutions
            j check1

        c2:
            feq.s t1, fa1, f0
            bnez t1, infinite_sols               # if true, infinite solutions
            beqz t1, no_sols                    #if true, no solutions
            j check2

        c3:
            feq.s t1, fa2, f0
            bnez t1, infinite_sols               # if true, infinite solutions
            beqz t1, no_sols                    #if true, no solutions
            j check3

        c4:
            feq.s t1, fa3, f0
            bnez t1, infinite_sols               # if true, infinite solutions
            beqz t1, no_sols                    #if true, no solutions
            j check4

        c5:
            feq.s t1, fa4, f0
            bnez t1, infinite_sols               # if true, infinite solutions
            beqz t1, no_sols                    #if true, no solutions
            j check5

    infinite_sols:
        li a0, 4
        la a1, infinite_string
        ecall

        j end_program
    
    no_sols:
        li a0, 4
        la a1, nosol_string
        ecall

        j end_program

    end_program:
        # Exit the program
        li a7, 10
        ecall