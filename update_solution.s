	.section	__TEXT,__text,regular,pure_instructions
	.macosx_version_min 10, 11
	.section	__TEXT,__literal8,8byte_literals
	.align	3
LCPI0_0:
	.quad	4607182418800017408     ## double 1
LCPI0_1:
	.quad	4696837146684686336     ## double 1.0E+6
LCPI0_2:
	.quad	4442235333156365461     ## double 9.9999999999999993E-12
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_newton_solve
	.align	4, 0x90
_newton_solve:                          ## @newton_solve
	.cfi_startproc
## BB#0:
	pushq	%rbp
Ltmp0:
	.cfi_def_cfa_offset 16
Ltmp1:
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
Ltmp2:
	.cfi_def_cfa_register %rbp
	pushq	%rbx
	subq	$7608, %rsp             ## imm = 0x1DB8
Ltmp3:
	.cfi_offset %rbx, -24
	movl	$2800, %eax             ## imm = 0xAF0
	movl	%eax, %r10d
	movq	___stack_chk_guard@GOTPCREL(%rip), %r11
	movq	(%r11), %r11
	movq	%r11, -16(%rbp)
	movq	%rdi, -7224(%rbp)
	movsd	%xmm0, -7232(%rbp)
	movq	%rsi, -7240(%rbp)
	movq	%rdx, -7248(%rbp)
	movq	%rcx, -7256(%rbp)
	movq	%r8, -7264(%rbp)
	movq	%r9, -7272(%rbp)
	movq	%r10, %rdi
	callq	_malloc
	movl	$2800, %ebx             ## imm = 0xAF0
	movl	%ebx, %edx
	movq	$-1, %rcx
	movq	%rax, -7280(%rbp)
	movq	-7280(%rbp), %rax
	movq	-7224(%rbp), %rsi
	movq	%rax, %rdi
	callq	___memcpy_chk
	leaq	L_.str(%rip), %rdi
	movq	%rax, -7360(%rbp)       ## 8-byte Spill
	movb	$0, %al
	callq	_printf
	leaq	L_.str.1(%rip), %rdi
	movsd	LCPI0_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movq	-7256(%rbp), %rcx
	movaps	%xmm0, %xmm1
	mulsd	(%rcx), %xmm1
	movq	-7256(%rbp), %rcx
	movaps	%xmm0, %xmm2
	mulsd	8(%rcx), %xmm2
	movq	-7256(%rbp), %rcx
	mulsd	16(%rcx), %xmm0
	movsd	%xmm0, -7368(%rbp)      ## 8-byte Spill
	movaps	%xmm1, %xmm0
	movaps	%xmm2, %xmm1
	movsd	-7368(%rbp), %xmm2      ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	movl	%eax, -7372(%rbp)       ## 4-byte Spill
	movb	$3, %al
	callq	_printf
	leaq	L_.str.2(%rip), %rdi
	movsd	LCPI0_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movq	-7256(%rbp), %rcx
	movaps	%xmm0, %xmm1
	mulsd	24(%rcx), %xmm1
	movq	-7256(%rbp), %rcx
	mulsd	32(%rcx), %xmm0
	movsd	%xmm0, -7384(%rbp)      ## 8-byte Spill
	movaps	%xmm1, %xmm0
	movsd	-7384(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	movl	%eax, -7388(%rbp)       ## 4-byte Spill
	movb	$2, %al
	callq	_printf
	leaq	L_.str.1(%rip), %rdi
	movq	-7256(%rbp), %rcx
	movsd	40(%rcx), %xmm0         ## xmm0 = mem[0],zero
	movq	-7256(%rbp), %rcx
	movsd	48(%rcx), %xmm1         ## xmm1 = mem[0],zero
	movq	-7256(%rbp), %rcx
	movsd	56(%rcx), %xmm2         ## xmm2 = mem[0],zero
	movl	%eax, -7392(%rbp)       ## 4-byte Spill
	movb	$3, %al
	callq	_printf
	leaq	L_.str.1(%rip), %rdi
	movq	-7256(%rbp), %rcx
	movsd	64(%rcx), %xmm0         ## xmm0 = mem[0],zero
	movq	-7256(%rbp), %rcx
	movsd	72(%rcx), %xmm1         ## xmm1 = mem[0],zero
	movq	-7256(%rbp), %rcx
	movsd	80(%rcx), %xmm2         ## xmm2 = mem[0],zero
	movl	%eax, -7396(%rbp)       ## 4-byte Spill
	movb	$3, %al
	callq	_printf
	leaq	L_.str.1(%rip), %rdi
	movq	-7256(%rbp), %rcx
	movsd	88(%rcx), %xmm0         ## xmm0 = mem[0],zero
	movq	-7256(%rbp), %rcx
	movsd	96(%rcx), %xmm1         ## xmm1 = mem[0],zero
	movq	-7256(%rbp), %rcx
	movsd	104(%rcx), %xmm2        ## xmm2 = mem[0],zero
	movl	%eax, -7400(%rbp)       ## 4-byte Spill
	movb	$3, %al
	callq	_printf
	leaq	L_.str.3(%rip), %rdi
	movsd	LCPI0_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movq	-7256(%rbp), %rcx
	movl	112(%rcx), %esi
	movq	-7256(%rbp), %rcx
	movaps	%xmm0, %xmm1
	mulsd	120(%rcx), %xmm1
	movq	-7256(%rbp), %rcx
	mulsd	128(%rcx), %xmm0
	movsd	%xmm0, -7408(%rbp)      ## 8-byte Spill
	movaps	%xmm1, %xmm0
	movsd	-7408(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	movl	%eax, -7412(%rbp)       ## 4-byte Spill
	movb	$2, %al
	callq	_printf
	movsd	LCPI0_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	leaq	-3616(%rbp), %rdi
	movq	-7224(%rbp), %rcx
	addq	$2400, %rcx             ## imm = 0x960
	movq	%rcx, %rsi
	movl	%eax, -7416(%rbp)       ## 4-byte Spill
	callq	_diff_coef
	movl	$0, -7284(%rbp)
LBB0_1:                                 ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB0_3 Depth 2
	cmpl	$3, -7284(%rbp)
	jge	LBB0_8
## BB#2:                                ##   in Loop: Header=BB0_1 Depth=1
	movl	$0, -7288(%rbp)
LBB0_3:                                 ##   Parent Loop BB0_1 Depth=1
                                        ## =>  This Inner Loop Header: Depth=2
	cmpl	$3, -7288(%rbp)
	jge	LBB0_6
## BB#4:                                ##   in Loop: Header=BB0_3 Depth=2
	leaq	L_.str.4(%rip), %rdi
	movl	-7284(%rbp), %esi
	movl	-7288(%rbp), %edx
	movb	$0, %al
	callq	_printf
	xorl	%edx, %edx
	movl	-7288(%rbp), %esi
	movl	-7284(%rbp), %ecx
	movl	%edx, %edi
	movl	%esi, -7420(%rbp)       ## 4-byte Spill
	movl	%edx, %esi
	movl	-7420(%rbp), %edx       ## 4-byte Reload
	movl	%eax, -7424(%rbp)       ## 4-byte Spill
	callq	_c_index
	xorl	%edi, %edi
	movl	$4, %esi
	movsd	LCPI0_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	shll	$1, %eax
	movslq	%eax, %r8
	mulsd	-3616(%rbp,%r8,8), %xmm0
	movl	-7288(%rbp), %edx
	movl	-7284(%rbp), %ecx
	movsd	%xmm0, -7432(%rbp)      ## 8-byte Spill
	callq	_c_index
	leaq	L_.str.5(%rip), %rdi
	movsd	LCPI0_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	mulsd	-3616(%rbp,%r8,8), %xmm0
	movsd	-7432(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	movsd	%xmm0, -7440(%rbp)      ## 8-byte Spill
	movaps	%xmm1, %xmm0
	movsd	-7440(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	movb	$2, %al
	callq	_printf
	movl	%eax, -7444(%rbp)       ## 4-byte Spill
## BB#5:                                ##   in Loop: Header=BB0_3 Depth=2
	movl	-7288(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7288(%rbp)
	jmp	LBB0_3
LBB0_6:                                 ##   in Loop: Header=BB0_1 Depth=1
	jmp	LBB0_7
LBB0_7:                                 ##   in Loop: Header=BB0_1 Depth=1
	movl	-7284(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7284(%rbp)
	jmp	LBB0_1
LBB0_8:
	movsd	LCPI0_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	leaq	-7216(%rbp), %rdi
	movq	-7224(%rbp), %rax
	addq	$2400, %rax             ## imm = 0x960
	movq	%rax, %rsi
	callq	_diff_coef
	movl	$0, -7292(%rbp)
LBB0_9:                                 ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB0_11 Depth 2
	cmpl	$3, -7292(%rbp)
	jge	LBB0_16
## BB#10:                               ##   in Loop: Header=BB0_9 Depth=1
	movl	$0, -7296(%rbp)
LBB0_11:                                ##   Parent Loop BB0_9 Depth=1
                                        ## =>  This Inner Loop Header: Depth=2
	cmpl	$3, -7296(%rbp)
	jge	LBB0_14
## BB#12:                               ##   in Loop: Header=BB0_11 Depth=2
	leaq	L_.str.6(%rip), %rdi
	movl	-7292(%rbp), %esi
	movl	-7296(%rbp), %edx
	movb	$0, %al
	callq	_printf
	xorl	%edx, %edx
	movl	-7296(%rbp), %esi
	movl	-7292(%rbp), %ecx
	movl	%edx, %edi
	movl	%esi, -7448(%rbp)       ## 4-byte Spill
	movl	%edx, %esi
	movl	-7448(%rbp), %edx       ## 4-byte Reload
	movl	%eax, -7452(%rbp)       ## 4-byte Spill
	callq	_c_index
	xorl	%edi, %edi
	movl	$4, %esi
	movsd	LCPI0_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	shll	$1, %eax
	movslq	%eax, %r8
	mulsd	-7216(%rbp,%r8,8), %xmm0
	movl	-7296(%rbp), %edx
	movl	-7292(%rbp), %ecx
	movsd	%xmm0, -7464(%rbp)      ## 8-byte Spill
	callq	_c_index
	leaq	L_.str.7(%rip), %rdi
	movsd	LCPI0_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	mulsd	-7216(%rbp,%r8,8), %xmm0
	movsd	-7464(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	movsd	%xmm0, -7472(%rbp)      ## 8-byte Spill
	movaps	%xmm1, %xmm0
	movsd	-7472(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	movb	$2, %al
	callq	_printf
	movl	%eax, -7476(%rbp)       ## 4-byte Spill
## BB#13:                               ##   in Loop: Header=BB0_11 Depth=2
	movl	-7296(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7296(%rbp)
	jmp	LBB0_11
LBB0_14:                                ##   in Loop: Header=BB0_9 Depth=1
	jmp	LBB0_15
LBB0_15:                                ##   in Loop: Header=BB0_9 Depth=1
	movl	-7292(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7292(%rbp)
	jmp	LBB0_9
LBB0_16:
	movl	$0, -7300(%rbp)
LBB0_17:                                ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB0_19 Depth 2
	cmpl	$3, -7300(%rbp)
	jge	LBB0_24
## BB#18:                               ##   in Loop: Header=BB0_17 Depth=1
	movl	$0, -7304(%rbp)
LBB0_19:                                ##   Parent Loop BB0_17 Depth=1
                                        ## =>  This Inner Loop Header: Depth=2
	cmpl	$3, -7304(%rbp)
	jge	LBB0_22
## BB#20:                               ##   in Loop: Header=BB0_19 Depth=2
	xorl	%eax, %eax
	movl	-7300(%rbp), %esi
	movl	-7304(%rbp), %edx
	movl	-7304(%rbp), %ecx
	movl	-7300(%rbp), %edi
	movl	%edi, -7480(%rbp)       ## 4-byte Spill
	movl	%eax, %edi
	movl	%esi, -7484(%rbp)       ## 4-byte Spill
	movl	%eax, %esi
	movl	%edx, -7488(%rbp)       ## 4-byte Spill
	movl	%ecx, %edx
	movl	-7480(%rbp), %ecx       ## 4-byte Reload
	callq	_c_index
	leaq	L_.str.8(%rip), %rdi
	movslq	%eax, %r8
	movq	-7224(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-7484(%rbp), %esi       ## 4-byte Reload
	movl	-7488(%rbp), %edx       ## 4-byte Reload
	movb	$1, %al
	callq	_printf
	movl	%eax, -7492(%rbp)       ## 4-byte Spill
## BB#21:                               ##   in Loop: Header=BB0_19 Depth=2
	movl	-7304(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7304(%rbp)
	jmp	LBB0_19
LBB0_22:                                ##   in Loop: Header=BB0_17 Depth=1
	jmp	LBB0_23
LBB0_23:                                ##   in Loop: Header=BB0_17 Depth=1
	movl	-7300(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7300(%rbp)
	jmp	LBB0_17
LBB0_24:
	movl	$0, -7308(%rbp)
LBB0_25:                                ## =>This Inner Loop Header: Depth=1
	cmpl	$3, -7308(%rbp)
	jge	LBB0_28
## BB#26:                               ##   in Loop: Header=BB0_25 Depth=1
	xorl	%eax, %eax
	movl	-7308(%rbp), %esi
	movl	-7308(%rbp), %edx
	movl	%eax, %edi
	movl	%esi, -7496(%rbp)       ## 4-byte Spill
	movl	%eax, %esi
	callq	_phi_index
	leaq	L_.str.9(%rip), %rdi
	movslq	%eax, %rcx
	movq	-7224(%rbp), %r8
	movsd	1800(%r8,%rcx,8), %xmm0 ## xmm0 = mem[0],zero
	movl	-7496(%rbp), %esi       ## 4-byte Reload
	movb	$1, %al
	callq	_printf
	movl	%eax, -7500(%rbp)       ## 4-byte Spill
## BB#27:                               ##   in Loop: Header=BB0_25 Depth=1
	movl	-7308(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7308(%rbp)
	jmp	LBB0_25
LBB0_28:
	movl	$225, %eax
	movl	%eax, %esi
	movq	-7224(%rbp), %rdi
	callq	_array_max
	movsd	LCPI0_0(%rip), %xmm1    ## xmm1 = mem[0],zero
	movsd	LCPI0_2(%rip), %xmm2    ## xmm2 = mem[0],zero
	mulsd	%xmm0, %xmm2
	movsd	%xmm2, -7320(%rbp)
	addsd	-7320(%rbp), %xmm1
	movsd	%xmm1, -7328(%rbp)
	movl	$0, -7332(%rbp)
	movl	$0, -7336(%rbp)
	movl	$0, -7340(%rbp)
## BB#29:
	cmpl	$10, -7340(%rbp)
	jge	LBB0_43
## BB#30:
	leaq	L_.str.10(%rip), %rdi
	movb	$0, %al
	callq	_printf
	leaq	L_.str.11(%rip), %rdi
	movsd	LCPI0_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movq	-7240(%rbp), %rcx
	movsd	(%rcx), %xmm1           ## xmm1 = mem[0],zero
	movq	-7240(%rbp), %rcx
	movsd	200(%rcx), %xmm2        ## xmm2 = mem[0],zero
	movq	-7240(%rbp), %rcx
	mulsd	400(%rcx), %xmm0
	movsd	%xmm0, -7512(%rbp)      ## 8-byte Spill
	movaps	%xmm1, %xmm0
	movaps	%xmm2, %xmm1
	movsd	-7512(%rbp), %xmm2      ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	movl	%eax, -7516(%rbp)       ## 4-byte Spill
	movb	$3, %al
	callq	_printf
	leaq	L_.str.12(%rip), %rdi
	movq	-7240(%rbp), %rcx
	movsd	600(%rcx), %xmm0        ## xmm0 = mem[0],zero
	movq	-7240(%rbp), %rcx
	movsd	800(%rcx), %xmm1        ## xmm1 = mem[0],zero
	movq	-7240(%rbp), %rcx
	movsd	1000(%rcx), %xmm2       ## xmm2 = mem[0],zero
	movl	%eax, -7520(%rbp)       ## 4-byte Spill
	movb	$3, %al
	callq	_printf
	leaq	L_.str.13(%rip), %rdi
	movq	-7240(%rbp), %rcx
	movsd	1200(%rcx), %xmm0       ## xmm0 = mem[0],zero
	movq	-7240(%rbp), %rcx
	movsd	1400(%rcx), %xmm1       ## xmm1 = mem[0],zero
	movl	%eax, -7524(%rbp)       ## 4-byte Spill
	movb	$2, %al
	callq	_printf
	leaq	L_.str.14(%rip), %rdi
	movq	-7240(%rbp), %rcx
	movsd	1600(%rcx), %xmm0       ## xmm0 = mem[0],zero
	movq	-7240(%rbp), %rcx
	movsd	1800(%rcx), %xmm1       ## xmm1 = mem[0],zero
	movq	-7240(%rbp), %rcx
	movsd	2000(%rcx), %xmm2       ## xmm2 = mem[0],zero
	movl	%eax, -7528(%rbp)       ## 4-byte Spill
	movb	$3, %al
	callq	_printf
	movq	-7272(%rbp), %rdi
	movq	-7224(%rbp), %rsi
	movq	-7280(%rbp), %rdx
	movq	-7240(%rbp), %rcx
	movq	-7248(%rbp), %r8
	movq	-7256(%rbp), %r9
	movl	%eax, -7532(%rbp)       ## 4-byte Spill
	callq	_ionmflux
	movl	$0, -7344(%rbp)
LBB0_31:                                ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB0_33 Depth 2
	cmpl	$3, -7344(%rbp)
	jge	LBB0_38
## BB#32:                               ##   in Loop: Header=BB0_31 Depth=1
	movl	$0, -7348(%rbp)
LBB0_33:                                ##   Parent Loop BB0_31 Depth=1
                                        ## =>  This Inner Loop Header: Depth=2
	cmpl	$3, -7348(%rbp)
	jge	LBB0_36
## BB#34:                               ##   in Loop: Header=BB0_33 Depth=2
	leaq	L_.str.15(%rip), %rdi
	movl	-7344(%rbp), %esi
	movl	-7348(%rbp), %edx
	movb	$0, %al
	callq	_printf
	xorl	%edx, %edx
	movl	-7348(%rbp), %esi
	movl	-7344(%rbp), %ecx
	movl	%edx, %edi
	movl	%esi, -7536(%rbp)       ## 4-byte Spill
	movl	%edx, %esi
	movl	-7536(%rbp), %edx       ## 4-byte Reload
	movl	%eax, -7540(%rbp)       ## 4-byte Spill
	callq	_c_index
	xorl	%ecx, %ecx
	movsd	LCPI0_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-7272(%rbp), %r9
	mulsd	(%r9,%r8,8), %xmm0
	movl	-7348(%rbp), %edx
	movl	-7344(%rbp), %eax
	movl	%ecx, %edi
	movl	%ecx, %esi
	movl	%eax, %ecx
	movsd	%xmm0, -7552(%rbp)      ## 8-byte Spill
	callq	_c_index
	xorl	%ecx, %ecx
	movslq	%eax, %r8
	movq	-7272(%rbp), %r9
	movsd	1800(%r9,%r8,8), %xmm1  ## xmm1 = mem[0],zero
	movl	-7348(%rbp), %edx
	movl	-7344(%rbp), %eax
	movl	%ecx, %edi
	movl	%ecx, %esi
	movl	%eax, %ecx
	movsd	%xmm1, -7560(%rbp)      ## 8-byte Spill
	callq	_c_index
	xorl	%ecx, %ecx
	movslq	%eax, %r8
	movq	-7272(%rbp), %r9
	movsd	3600(%r9,%r8,8), %xmm2  ## xmm2 = mem[0],zero
	movl	-7348(%rbp), %edx
	movl	-7344(%rbp), %eax
	movl	%ecx, %edi
	movl	%ecx, %esi
	movl	%eax, %ecx
	movsd	%xmm2, -7568(%rbp)      ## 8-byte Spill
	callq	_c_index
	leaq	L_.str.16(%rip), %rdi
	movslq	%eax, %r8
	movq	-7272(%rbp), %r9
	movsd	5400(%r9,%r8,8), %xmm3  ## xmm3 = mem[0],zero
	movsd	-7552(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	movsd	-7560(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	movsd	-7568(%rbp), %xmm2      ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	movb	$4, %al
	callq	_printf
	movl	%eax, -7572(%rbp)       ## 4-byte Spill
## BB#35:                               ##   in Loop: Header=BB0_33 Depth=2
	movl	-7348(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7348(%rbp)
	jmp	LBB0_33
LBB0_36:                                ##   in Loop: Header=BB0_31 Depth=1
	jmp	LBB0_37
LBB0_37:                                ##   in Loop: Header=BB0_31 Depth=1
	movl	-7344(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7344(%rbp)
	jmp	LBB0_31
LBB0_38:
	movq	-7272(%rbp), %rdi
	movq	-7224(%rbp), %rsi
	movq	-7256(%rbp), %rdx
	callq	_wflowm
	movl	$0, -7352(%rbp)
LBB0_39:                                ## =>This Inner Loop Header: Depth=1
	cmpl	$2, -7352(%rbp)
	jge	LBB0_42
## BB#40:                               ##   in Loop: Header=BB0_39 Depth=1
	leaq	L_.str.17(%rip), %rdi
	movl	-7352(%rbp), %esi
	movb	$0, %al
	callq	_printf
	movl	-7332(%rbp), %edi
	movl	-7336(%rbp), %esi
	movl	-7352(%rbp), %edx
	movl	%eax, -7576(%rbp)       ## 4-byte Spill
	callq	_al_index
	movslq	%eax, %rcx
	movq	-7272(%rbp), %r8
	movsd	7200(%r8,%rcx,8), %xmm0 ## xmm0 = mem[0],zero
	movl	-7332(%rbp), %edi
	movl	-7336(%rbp), %esi
	movl	-7352(%rbp), %edx
	movsd	%xmm0, -7584(%rbp)      ## 8-byte Spill
	callq	_al_index
	movslq	%eax, %rcx
	movq	-7272(%rbp), %r8
	movsd	7600(%r8,%rcx,8), %xmm1 ## xmm1 = mem[0],zero
	movl	-7332(%rbp), %edi
	movl	-7336(%rbp), %esi
	movl	-7352(%rbp), %edx
	movsd	%xmm1, -7592(%rbp)      ## 8-byte Spill
	callq	_al_index
	leaq	L_.str.18(%rip), %rdi
	movslq	%eax, %rcx
	movq	-7272(%rbp), %r8
	movsd	8000(%r8,%rcx,8), %xmm2 ## xmm2 = mem[0],zero
	movsd	-7584(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	movsd	-7592(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	movb	$3, %al
	callq	_printf
	movl	%eax, -7596(%rbp)       ## 4-byte Spill
## BB#41:                               ##   in Loop: Header=BB0_39 Depth=1
	movl	-7352(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -7352(%rbp)
	jmp	LBB0_39
LBB0_42:
	leaq	-7216(%rbp), %r8
	leaq	-3616(%rbp), %rcx
	movq	-7264(%rbp), %rax
	movq	29528(%rax), %rdi
	movq	-7280(%rbp), %rsi
	movq	-7224(%rbp), %rdx
	movsd	-7232(%rbp), %xmm0      ## xmm0 = mem[0],zero
	movq	-7272(%rbp), %r9
	movq	-7256(%rbp), %rax
	movq	%rax, (%rsp)
	callq	_calc_residual
	leaq	-7216(%rbp), %r8
	leaq	-3616(%rbp), %rcx
	movq	-7264(%rbp), %rdx
	movq	29536(%rdx), %rdi
	movq	-7280(%rbp), %rsi
	movq	-7224(%rbp), %rdx
	movsd	-7232(%rbp), %xmm0      ## xmm0 = mem[0],zero
	movq	-7272(%rbp), %r9
	movq	-7256(%rbp), %r10
	movq	%r10, (%rsp)
	movl	%eax, -7600(%rbp)       ## 4-byte Spill
	callq	_calc_jacobian
	movl	%eax, -7604(%rbp)       ## 4-byte Spill
LBB0_43:
	movsd	-7328(%rbp), %xmm0      ## xmm0 = mem[0],zero
	ucomisd	-7320(%rbp), %xmm0
	jbe	LBB0_45
## BB#44:
	leaq	L_.str.19(%rip), %rsi
	movq	___stderrp@GOTPCREL(%rip), %rax
	movq	(%rax), %rdi
	movb	$0, %al
	callq	_fprintf
	movl	$1, %edi
	movl	%eax, -7608(%rbp)       ## 4-byte Spill
	callq	_exit
LBB0_45:
	movq	___stack_chk_guard@GOTPCREL(%rip), %rax
	movq	(%rax), %rax
	cmpq	-16(%rbp), %rax
	jne	LBB0_47
## BB#46:
	addq	$7608, %rsp             ## imm = 0x1DB8
	popq	%rbx
	popq	%rbp
	retq
LBB0_47:
	callq	___stack_chk_fail
	.cfi_endproc

	.section	__TEXT,__literal8,8byte_literals
	.align	3
LCPI1_0:
	.quad	4607182418800017408     ## double 1
LCPI1_1:
	.quad	4576918229304087675     ## double 0.01
LCPI1_2:
	.quad	4611686018427387904     ## double 2
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_calc_residual
	.align	4, 0x90
_calc_residual:                         ## @calc_residual
	.cfi_startproc
## BB#0:
	pushq	%rbp
Ltmp4:
	.cfi_def_cfa_offset 16
Ltmp5:
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
Ltmp6:
	.cfi_def_cfa_register %rbp
	subq	$1072, %rsp             ## imm = 0x430
	movq	16(%rbp), %rax
	movq	%rdi, -16(%rbp)
	movq	%rsi, -24(%rbp)
	movq	%rdx, -32(%rbp)
	movsd	%xmm0, -40(%rbp)
	movq	%rcx, -48(%rbp)
	movq	%r8, -56(%rbp)
	movq	%r9, -64(%rbp)
	movq	%rax, -72(%rbp)
	movq	-32(%rbp), %rax
	movq	%rax, -80(%rbp)
	movq	-32(%rbp), %rax
	addq	$1800, %rax             ## imm = 0x708
	movq	%rax, -88(%rbp)
	movq	-32(%rbp), %rax
	addq	$2400, %rax             ## imm = 0x960
	movq	%rax, -96(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, -104(%rbp)
	movq	-24(%rbp), %rax
	addq	$2400, %rax             ## imm = 0x960
	movq	%rax, -112(%rbp)
	movl	$0, -180(%rbp)
LBB1_1:                                 ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB1_3 Depth 2
                                        ##       Child Loop BB1_5 Depth 3
                                        ##         Child Loop BB1_7 Depth 4
	cmpl	$5, -180(%rbp)
	jge	LBB1_44
## BB#2:                                ##   in Loop: Header=BB1_1 Depth=1
	movl	$0, -184(%rbp)
LBB1_3:                                 ##   Parent Loop BB1_1 Depth=1
                                        ## =>  This Loop Header: Depth=2
                                        ##       Child Loop BB1_5 Depth 3
                                        ##         Child Loop BB1_7 Depth 4
	cmpl	$5, -184(%rbp)
	jge	LBB1_42
## BB#4:                                ##   in Loop: Header=BB1_3 Depth=2
	movl	$0, -172(%rbp)
LBB1_5:                                 ##   Parent Loop BB1_1 Depth=1
                                        ##     Parent Loop BB1_3 Depth=2
                                        ## =>    This Loop Header: Depth=3
                                        ##         Child Loop BB1_7 Depth 4
	cmpl	$3, -172(%rbp)
	jge	LBB1_40
## BB#6:                                ##   in Loop: Header=BB1_5 Depth=3
	movl	$0, -176(%rbp)
LBB1_7:                                 ##   Parent Loop BB1_1 Depth=1
                                        ##     Parent Loop BB1_3 Depth=2
                                        ##       Parent Loop BB1_5 Depth=3
                                        ## =>      This Inner Loop Header: Depth=4
	cmpl	$2, -176(%rbp)
	jge	LBB1_22
## BB#8:                                ##   in Loop: Header=BB1_7 Depth=4
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -120(%rbp)
	movsd	%xmm0, -144(%rbp)
	cmpl	$0, -180(%rbp)
	jle	LBB1_10
## BB#9:                                ##   in Loop: Header=BB1_7 Depth=4
	movl	-180(%rbp), %eax
	subl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %eax
	subl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -200(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -208(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	-208(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm1
	movsd	-200(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -120(%rbp)
	movsd	-120(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -216(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	movl	-180(%rbp), %eax
	subl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -224(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	leaq	_z(%rip), %r8
	movsd	-224(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm1, -232(%rbp)       ## 8-byte Spill
	movsd	%xmm0, -240(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %eax
	subl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	%eax, %edi
	movsd	%xmm0, -248(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movsd	LCPI1_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-248(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm1
	movsd	-240(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	movsd	-232(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm2, %xmm1
	movsd	-216(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	mulsd	-40(%rbp), %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -120(%rbp)
LBB1_10:                                ##   in Loop: Header=BB1_7 Depth=4
	cmpl	$4, -180(%rbp)
	jge	LBB1_12
## BB#11:                               ##   in Loop: Header=BB1_7 Depth=4
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -256(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %eax
	addl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -264(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	-264(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm1
	movsd	-256(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -144(%rbp)
	movsd	-144(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-180(%rbp), %eax
	addl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -272(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -280(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	leaq	_z(%rip), %r8
	movsd	-280(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-180(%rbp), %eax
	addl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	%eax, %edi
	movsd	%xmm1, -288(%rbp)       ## 8-byte Spill
	movsd	%xmm0, -296(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm0, -304(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movsd	LCPI1_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-304(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm1
	movsd	-296(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	movsd	-288(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm2, %xmm1
	movsd	-272(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	mulsd	-40(%rbp), %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -144(%rbp)
LBB1_12:                                ##   in Loop: Header=BB1_7 Depth=4
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -128(%rbp)
	movsd	%xmm0, -152(%rbp)
	cmpl	$0, -184(%rbp)
	jle	LBB1_14
## BB#13:                               ##   in Loop: Header=BB1_7 Depth=4
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	subl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	callq	_c_index
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	subl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -312(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -320(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	-320(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm1
	movsd	-312(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -128(%rbp)
	movsd	-128(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -328(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	subl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -336(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	leaq	_z(%rip), %r8
	movsd	-336(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm1, -344(%rbp)       ## 8-byte Spill
	movsd	%xmm0, -352(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	subl	$1, %eax
	movl	-176(%rbp), %edx
	movl	%eax, %esi
	movsd	%xmm0, -360(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movsd	LCPI1_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-360(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm1
	movsd	-352(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	movsd	-344(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm2, %xmm1
	movsd	-328(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	mulsd	-40(%rbp), %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -128(%rbp)
LBB1_14:                                ##   in Loop: Header=BB1_7 Depth=4
	cmpl	$4, -184(%rbp)
	jge	LBB1_16
## BB#15:                               ##   in Loop: Header=BB1_7 Depth=4
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -368(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	addl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -376(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	-376(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm1
	movsd	-368(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -152(%rbp)
	movsd	-152(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	addl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -384(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -392(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	leaq	_z(%rip), %r8
	movsd	-392(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	addl	$1, %eax
	movl	-176(%rbp), %edx
	movl	%eax, %esi
	movsd	%xmm1, -400(%rbp)       ## 8-byte Spill
	movsd	%xmm0, -408(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm0, -416(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movsd	LCPI1_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-416(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm1
	movsd	-408(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	movsd	-400(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm2, %xmm1
	movsd	-384(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	mulsd	-40(%rbp), %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -152(%rbp)
LBB1_16:                                ##   in Loop: Header=BB1_7 Depth=4
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	callq	_al_index
	movslq	%eax, %rcx
	movq	-96(%rbp), %r8
	movsd	(%r8,%rcx,8), %xmm0     ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -424(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	-424(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	(%r9,%r8,8), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm0, -432(%rbp)       ## 8-byte Spill
	callq	_al_index
	movslq	%eax, %r8
	movq	-112(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -440(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	-440(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	(%r9,%r8,8), %xmm0
	movsd	-432(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movsd	%xmm1, -136(%rbp)
	movsd	-120(%rbp), %xmm0       ## xmm0 = mem[0],zero
	subsd	-144(%rbp), %xmm0
	addsd	-128(%rbp), %xmm0
	subsd	-152(%rbp), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -448(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	movsd	-448(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm0, %xmm1
	addsd	-136(%rbp), %xmm1
	movsd	%xmm1, -136(%rbp)
	movq	-16(%rbp), %rdi
	movl	-180(%rbp), %eax
	movl	-184(%rbp), %esi
	movl	-172(%rbp), %edx
	movl	-176(%rbp), %ecx
	movq	%rdi, -456(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$1, %edx
	movsd	-136(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-456(%rbp), %rdi        ## 8-byte Reload
	movl	%eax, %esi
	callq	_VecSetValue
	movl	%eax, -188(%rbp)
## BB#17:                               ##   in Loop: Header=BB1_7 Depth=4
	cmpl	$0, -188(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB1_19
## BB#18:
	movl	$1, %eax
	movl	$231, %esi
	leaq	L___func__.calc_residual(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-188(%rbp), %r8d
	movq	%rdi, -464(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-464(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB1_97
LBB1_19:                                ##   in Loop: Header=BB1_7 Depth=4
	jmp	LBB1_20
LBB1_20:                                ##   in Loop: Header=BB1_7 Depth=4
	jmp	LBB1_21
LBB1_21:                                ##   in Loop: Header=BB1_7 Depth=4
	movl	-176(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -176(%rbp)
	jmp	LBB1_7
LBB1_22:                                ##   in Loop: Header=BB1_5 Depth=3
	xorl	%edx, %edx
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	callq	_al_index
	movl	$1, %edx
	movsd	LCPI1_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %rcx
	movq	-96(%rbp), %r8
	subsd	(%r8,%rcx,8), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movsd	%xmm0, -472(%rbp)       ## 8-byte Spill
	callq	_al_index
	xorl	%edx, %edx
	movslq	%eax, %rcx
	movq	-96(%rbp), %r8
	movsd	-472(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r8,%rcx,8), %xmm0
	movsd	%xmm0, -160(%rbp)
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	callq	_al_index
	movl	$1, %edx
	movsd	LCPI1_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %rcx
	movq	-96(%rbp), %r8
	subsd	(%r8,%rcx,8), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movsd	%xmm0, -480(%rbp)       ## 8-byte Spill
	callq	_al_index
	xorps	%xmm0, %xmm0
	movslq	%eax, %rcx
	movq	-96(%rbp), %r8
	movsd	-480(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	(%r8,%rcx,8), %xmm1
	movsd	%xmm1, -168(%rbp)
	movl	$2, -176(%rbp)
	movsd	%xmm0, -120(%rbp)
	movsd	%xmm0, -144(%rbp)
	cmpl	$0, -180(%rbp)
	jle	LBB1_24
## BB#23:                               ##   in Loop: Header=BB1_5 Depth=3
	movl	-180(%rbp), %eax
	subl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %eax
	subl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -488(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -496(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	-496(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm1
	movsd	-488(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -120(%rbp)
	movsd	-120(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -504(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	movl	-180(%rbp), %eax
	subl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -512(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	leaq	_z(%rip), %r8
	movsd	-512(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm1, -520(%rbp)       ## 8-byte Spill
	movsd	%xmm0, -528(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %eax
	subl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	%eax, %edi
	movsd	%xmm0, -536(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movsd	LCPI1_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-536(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm1
	movsd	-528(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	movsd	-520(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm2, %xmm1
	movsd	-504(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	mulsd	-40(%rbp), %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -120(%rbp)
LBB1_24:                                ##   in Loop: Header=BB1_5 Depth=3
	cmpl	$4, -180(%rbp)
	jge	LBB1_26
## BB#25:                               ##   in Loop: Header=BB1_5 Depth=3
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -544(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %eax
	addl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -552(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	-552(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm1
	movsd	-544(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -144(%rbp)
	movsd	-144(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-180(%rbp), %eax
	addl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -560(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -568(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	leaq	_z(%rip), %r8
	movsd	-568(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-180(%rbp), %eax
	addl	$1, %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	%eax, %edi
	movsd	%xmm1, -576(%rbp)       ## 8-byte Spill
	movsd	%xmm0, -584(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm0, -592(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movsd	LCPI1_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-592(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm1
	movsd	-584(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	movsd	-576(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm2, %xmm1
	movsd	-560(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	mulsd	-40(%rbp), %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -144(%rbp)
LBB1_26:                                ##   in Loop: Header=BB1_5 Depth=3
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -128(%rbp)
	movsd	%xmm0, -152(%rbp)
	cmpl	$0, -184(%rbp)
	jle	LBB1_28
## BB#27:                               ##   in Loop: Header=BB1_5 Depth=3
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	subl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	callq	_c_index
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	subl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -600(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -608(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	-608(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm1
	movsd	-600(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -128(%rbp)
	movsd	-128(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -616(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	subl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -624(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	leaq	_z(%rip), %r8
	movsd	-624(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm1, -632(%rbp)       ## 8-byte Spill
	movsd	%xmm0, -640(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	subl	$1, %eax
	movl	-176(%rbp), %edx
	movl	%eax, %esi
	movsd	%xmm0, -648(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movsd	LCPI1_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-648(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm1
	movsd	-640(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	movsd	-632(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm2, %xmm1
	movsd	-616(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	mulsd	-40(%rbp), %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -128(%rbp)
LBB1_28:                                ##   in Loop: Header=BB1_5 Depth=3
	cmpl	$4, -184(%rbp)
	jge	LBB1_30
## BB#29:                               ##   in Loop: Header=BB1_5 Depth=3
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -656(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	addl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -664(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-104(%rbp), %r9
	movsd	-664(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm1
	movsd	-656(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -152(%rbp)
	movsd	-152(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	addl	$1, %eax
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -672(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -680(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	leaq	_z(%rip), %r8
	movsd	-680(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %eax
	addl	$1, %eax
	movl	-176(%rbp), %edx
	movl	%eax, %esi
	movsd	%xmm1, -688(%rbp)       ## 8-byte Spill
	movsd	%xmm0, -696(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm0, -704(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movsd	LCPI1_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-704(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm1
	movsd	-696(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	movsd	-688(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm2, %xmm1
	movsd	-672(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	mulsd	-40(%rbp), %xmm2
	divsd	%xmm0, %xmm2
	movsd	%xmm2, -152(%rbp)
LBB1_30:                                ##   in Loop: Header=BB1_5 Depth=3
	movsd	-160(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -712(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	-712(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	(%r9,%r8,8), %xmm0
	movsd	-168(%rbp), %xmm1       ## xmm1 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -720(%rbp)       ## 8-byte Spill
	movsd	%xmm1, -728(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	-728(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	(%r9,%r8,8), %xmm0
	movsd	-720(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movsd	%xmm1, -136(%rbp)
	movsd	-120(%rbp), %xmm0       ## xmm0 = mem[0],zero
	subsd	-144(%rbp), %xmm0
	addsd	-128(%rbp), %xmm0
	subsd	-152(%rbp), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -736(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	movsd	-736(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm0, %xmm1
	addsd	-136(%rbp), %xmm1
	movsd	%xmm1, -136(%rbp)
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-56(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movsd	LCPI1_2(%rip), %xmm1    ## xmm1 = mem[0],zero
	callq	_pow
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -744(%rbp)       ## 8-byte Spill
	callq	_c_index
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-56(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movsd	LCPI1_2(%rip), %xmm1    ## xmm1 = mem[0],zero
	callq	_pow
	movsd	-744(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm0, %xmm1
	sqrtsd	%xmm1, %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -752(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	leaq	_cbath(%rip), %r8
	movslq	%eax, %r9
	movq	-104(%rbp), %r10
	movsd	(%r10,%r9,8), %xmm1     ## xmm1 = mem[0],zero
	movslq	-172(%rbp), %r9
	addsd	(%r8,%r9,8), %xmm1
	movsd	-752(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm2, -760(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	callq	_log
	leaq	_cbath(%rip), %r8
	movslq	-172(%rbp), %r9
	movsd	(%r8,%r9,8), %xmm1      ## xmm1 = mem[0],zero
	movsd	%xmm0, -768(%rbp)       ## 8-byte Spill
	movaps	%xmm1, %xmm0
	callq	_log
	leaq	_z(%rip), %r8
	movsd	-768(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm1, -776(%rbp)       ## 8-byte Spill
	movsd	%xmm0, -784(%rbp)       ## 8-byte Spill
	callq	_phi_index
	xorps	%xmm0, %xmm0
	leaq	_z(%rip), %r8
	movslq	%eax, %r9
	movq	-88(%rbp), %r10
	movsd	-784(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	mulsd	(%r10,%r9,8), %xmm1
	movsd	-776(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	addsd	%xmm1, %xmm2
	movslq	-172(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm1
	mulsd	%xmm0, %xmm1
	subsd	%xmm1, %xmm2
	movsd	-760(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	%xmm2, %xmm0
	mulsd	-40(%rbp), %xmm0
	addsd	-136(%rbp), %xmm0
	movsd	%xmm0, -136(%rbp)
	cmpl	$0, -180(%rbp)
	sete	%r11b
	andb	$1, %r11b
	movzbl	%r11b, %eax
	cmpl	$4, -184(%rbp)
	sete	%r11b
	andb	$1, %r11b
	movzbl	%r11b, %ecx
	andl	%ecx, %eax
	cmpl	$0, -172(%rbp)
	sete	%r11b
	andb	$1, %r11b
	movzbl	%r11b, %ecx
	andl	%ecx, %eax
	cmpl	$0, %eax
	je	LBB1_32
## BB#31:                               ##   in Loop: Header=BB1_5 Depth=3
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	callq	_c_index
	leaq	L_.str.22(%rip), %rdi
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-56(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movb	$1, %al
	callq	_printf
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movl	%eax, -788(%rbp)        ## 4-byte Spill
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-56(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movsd	LCPI1_2(%rip), %xmm1    ## xmm1 = mem[0],zero
	callq	_pow
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	%xmm0, -800(%rbp)       ## 8-byte Spill
	callq	_c_index
	leaq	L_.str.22(%rip), %rdi
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-56(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movsd	LCPI1_2(%rip), %xmm1    ## xmm1 = mem[0],zero
	movq	%rdi, -808(%rbp)        ## 8-byte Spill
	callq	_pow
	movsd	-800(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm0, %xmm1
	sqrtsd	%xmm1, %xmm0
	movq	-808(%rbp), %rdi        ## 8-byte Reload
	movb	$1, %al
	callq	_printf
	leaq	L_.str.23(%rip), %rdi
	movl	-180(%rbp), %esi
	movl	-184(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	-136(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	%eax, -812(%rbp)        ## 4-byte Spill
	movb	$1, %al
	callq	_printf
	movl	%eax, -816(%rbp)        ## 4-byte Spill
LBB1_32:                                ##   in Loop: Header=BB1_5 Depth=3
	cmpl	$4, -180(%rbp)
	sete	%al
	andb	$1, %al
	movzbl	%al, %ecx
	cmpl	$3, -184(%rbp)
	sete	%al
	andb	$1, %al
	movzbl	%al, %edx
	andl	%edx, %ecx
	cmpl	$0, -172(%rbp)
	sete	%al
	andb	$1, %al
	movzbl	%al, %edx
	andl	%edx, %ecx
	cmpl	$0, %ecx
	je	LBB1_34
## BB#33:                               ##   in Loop: Header=BB1_5 Depth=3
	leaq	L_.str.23(%rip), %rdi
	movl	-180(%rbp), %esi
	movl	-184(%rbp), %edx
	movl	-172(%rbp), %ecx
	movsd	-136(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movb	$1, %al
	callq	_printf
	movl	%eax, -820(%rbp)        ## 4-byte Spill
LBB1_34:                                ##   in Loop: Header=BB1_5 Depth=3
	movq	-16(%rbp), %rdi
	movl	-180(%rbp), %eax
	movl	-184(%rbp), %esi
	movl	-172(%rbp), %edx
	movl	-176(%rbp), %ecx
	movq	%rdi, -832(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$1, %edx
	movsd	-136(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-832(%rbp), %rdi        ## 8-byte Reload
	movl	%eax, %esi
	callq	_VecSetValue
	movl	%eax, -188(%rbp)
## BB#35:                               ##   in Loop: Header=BB1_5 Depth=3
	cmpl	$0, -188(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB1_37
## BB#36:
	movl	$1, %eax
	movl	$281, %esi              ## imm = 0x119
	leaq	L___func__.calc_residual(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-188(%rbp), %r8d
	movq	%rdi, -840(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-840(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB1_97
LBB1_37:                                ##   in Loop: Header=BB1_5 Depth=3
	jmp	LBB1_38
LBB1_38:                                ##   in Loop: Header=BB1_5 Depth=3
	jmp	LBB1_39
LBB1_39:                                ##   in Loop: Header=BB1_5 Depth=3
	movl	-172(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -172(%rbp)
	jmp	LBB1_5
LBB1_40:                                ##   in Loop: Header=BB1_3 Depth=2
	jmp	LBB1_41
LBB1_41:                                ##   in Loop: Header=BB1_3 Depth=2
	movl	-184(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -184(%rbp)
	jmp	LBB1_3
LBB1_42:                                ##   in Loop: Header=BB1_1 Depth=1
	jmp	LBB1_43
LBB1_43:                                ##   in Loop: Header=BB1_1 Depth=1
	movl	-180(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -180(%rbp)
	jmp	LBB1_1
LBB1_44:
	movl	$0, -180(%rbp)
LBB1_45:                                ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB1_47 Depth 2
                                        ##       Child Loop BB1_49 Depth 3
                                        ##       Child Loop BB1_61 Depth 3
	cmpl	$5, -180(%rbp)
	jge	LBB1_68
## BB#46:                               ##   in Loop: Header=BB1_45 Depth=1
	movl	$0, -184(%rbp)
LBB1_47:                                ##   Parent Loop BB1_45 Depth=1
                                        ## =>  This Loop Header: Depth=2
                                        ##       Child Loop BB1_49 Depth 3
                                        ##       Child Loop BB1_61 Depth 3
	cmpl	$5, -184(%rbp)
	jge	LBB1_66
## BB#48:                               ##   in Loop: Header=BB1_47 Depth=2
	movl	$0, -176(%rbp)
LBB1_49:                                ##   Parent Loop BB1_45 Depth=1
                                        ##     Parent Loop BB1_47 Depth=2
                                        ## =>    This Inner Loop Header: Depth=3
	cmpl	$2, -176(%rbp)
	jge	LBB1_56
## BB#50:                               ##   in Loop: Header=BB1_49 Depth=3
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	callq	_al_index
	leaq	_z(%rip), %rsi
	movslq	%eax, %rcx
	movq	-96(%rbp), %r8
	movsd	(%r8,%rcx,8), %xmm0     ## xmm0 = mem[0],zero
	movq	-80(%rbp), %rdi
	movl	-180(%rbp), %edx
	movl	-184(%rbp), %ecx
	movl	-176(%rbp), %r8d
	movsd	%xmm0, -848(%rbp)       ## 8-byte Spill
	callq	_cz
	xorl	%eax, %eax
	movsd	-848(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	mulsd	%xmm0, %xmm1
	movl	-176(%rbp), %edx
	movl	%eax, %edi
	movl	%eax, %esi
	movsd	%xmm1, -856(%rbp)       ## 8-byte Spill
	callq	_phi_index
	xorl	%ecx, %ecx
	movslq	%eax, %r9
	movq	-72(%rbp), %r10
	movsd	64(%r10,%r9,8), %xmm0   ## xmm0 = mem[0],zero
	movl	-176(%rbp), %edx
	movl	%ecx, %edi
	movl	%ecx, %esi
	movsd	%xmm0, -864(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movl	$3, %edx
	movslq	%eax, %r9
	movq	-72(%rbp), %r10
	movsd	-864(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	40(%r10,%r9,8), %xmm0
	movsd	-856(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm0, %xmm1
	movsd	%xmm1, -136(%rbp)
	movq	-16(%rbp), %rdi
	movl	-180(%rbp), %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %ecx
	movq	%rdi, -872(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$1, %edx
	movsd	-136(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-872(%rbp), %rdi        ## 8-byte Reload
	movl	%eax, %esi
	callq	_VecSetValue
	movl	%eax, -188(%rbp)
## BB#51:                               ##   in Loop: Header=BB1_49 Depth=3
	cmpl	$0, -188(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB1_53
## BB#52:
	movl	$1, %eax
	movl	$296, %esi              ## imm = 0x128
	leaq	L___func__.calc_residual(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-188(%rbp), %r8d
	movq	%rdi, -880(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-880(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB1_97
LBB1_53:                                ##   in Loop: Header=BB1_49 Depth=3
	jmp	LBB1_54
LBB1_54:                                ##   in Loop: Header=BB1_49 Depth=3
	jmp	LBB1_55
LBB1_55:                                ##   in Loop: Header=BB1_49 Depth=3
	movl	-176(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -176(%rbp)
	jmp	LBB1_49
LBB1_56:                                ##   in Loop: Header=BB1_47 Depth=2
	xorl	%edx, %edx
	movl	$2, -176(%rbp)
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	callq	_al_index
	movl	$1, %edx
	movsd	LCPI1_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %rcx
	movq	-96(%rbp), %r8
	subsd	(%r8,%rcx,8), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movsd	%xmm0, -888(%rbp)       ## 8-byte Spill
	callq	_al_index
	leaq	_z(%rip), %rsi
	movslq	%eax, %rcx
	movq	-96(%rbp), %r8
	movsd	-888(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r8,%rcx,8), %xmm0
	movq	-80(%rbp), %rdi
	movl	-180(%rbp), %edx
	movl	-184(%rbp), %ecx
	movl	-176(%rbp), %r8d
	movsd	%xmm0, -896(%rbp)       ## 8-byte Spill
	callq	_cz
	xorl	%eax, %eax
	movsd	-896(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	mulsd	%xmm0, %xmm1
	movl	-176(%rbp), %edx
	movl	%eax, %edi
	movl	%eax, %esi
	movsd	%xmm1, -904(%rbp)       ## 8-byte Spill
	callq	_phi_index
	xorl	%ecx, %ecx
	movslq	%eax, %r9
	movq	-72(%rbp), %r10
	movsd	64(%r10,%r9,8), %xmm0   ## xmm0 = mem[0],zero
	movl	-176(%rbp), %edx
	movl	%ecx, %edi
	movl	%ecx, %esi
	movsd	%xmm0, -912(%rbp)       ## 8-byte Spill
	callq	_phi_index
	movl	$3, %edx
	movslq	%eax, %r9
	movq	-72(%rbp), %r10
	movsd	-912(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	40(%r10,%r9,8), %xmm0
	movsd	-904(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm0, %xmm1
	movsd	%xmm1, -136(%rbp)
	movq	-16(%rbp), %rdi
	movl	-180(%rbp), %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %ecx
	movq	%rdi, -920(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$1, %edx
	movsd	-136(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-920(%rbp), %rdi        ## 8-byte Reload
	movl	%eax, %esi
	callq	_VecSetValue
	movl	%eax, -188(%rbp)
## BB#57:                               ##   in Loop: Header=BB1_47 Depth=2
	cmpl	$0, -188(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB1_59
## BB#58:
	movl	$1, %eax
	movl	$301, %esi              ## imm = 0x12D
	leaq	L___func__.calc_residual(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-188(%rbp), %r8d
	movq	%rdi, -928(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-928(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB1_97
LBB1_59:                                ##   in Loop: Header=BB1_47 Depth=2
	jmp	LBB1_60
LBB1_60:                                ##   in Loop: Header=BB1_47 Depth=2
	movl	$0, -176(%rbp)
LBB1_61:                                ##   Parent Loop BB1_45 Depth=1
                                        ##     Parent Loop BB1_47 Depth=2
                                        ## =>    This Inner Loop Header: Depth=3
	cmpl	$2, -176(%rbp)
	jge	LBB1_64
## BB#62:                               ##   in Loop: Header=BB1_61 Depth=3
	movl	$4, %edx
	movq	-16(%rbp), %rdi
	movl	-180(%rbp), %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %ecx
	movq	%rdi, -936(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	%eax, -940(%rbp)        ## 4-byte Spill
	callq	_al_index
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm0, -952(%rbp)       ## 8-byte Spill
	callq	_al_index
	movslq	%eax, %r8
	movq	-112(%rbp), %r9
	movsd	-952(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm0, -960(%rbp)       ## 8-byte Spill
	callq	_al_index
	movl	$1, %edx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	7200(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	movsd	-960(%rbp), %xmm1       ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm0, %xmm1
	movq	-936(%rbp), %rdi        ## 8-byte Reload
	movl	-940(%rbp), %esi        ## 4-byte Reload
	movaps	%xmm1, %xmm0
	callq	_VecSetValue
	movl	%eax, -188(%rbp)
## BB#63:                               ##   in Loop: Header=BB1_61 Depth=3
	movl	-176(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -176(%rbp)
	jmp	LBB1_61
LBB1_64:                                ##   in Loop: Header=BB1_47 Depth=2
	jmp	LBB1_65
LBB1_65:                                ##   in Loop: Header=BB1_47 Depth=2
	movl	-184(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -184(%rbp)
	jmp	LBB1_47
LBB1_66:                                ##   in Loop: Header=BB1_45 Depth=1
	jmp	LBB1_67
LBB1_67:                                ##   in Loop: Header=BB1_45 Depth=1
	movl	-180(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -180(%rbp)
	jmp	LBB1_45
LBB1_68:
	movq	-16(%rbp), %rdi
	callq	_VecAssemblyBegin
	movl	%eax, -188(%rbp)
## BB#69:
	cmpl	$0, -188(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB1_71
## BB#70:
	movl	$1, %eax
	movl	$314, %esi              ## imm = 0x13A
	leaq	L___func__.calc_residual(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-188(%rbp), %r8d
	movq	%rdi, -968(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-968(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB1_97
LBB1_71:
	jmp	LBB1_72
LBB1_72:
	movq	-16(%rbp), %rdi
	callq	_VecAssemblyEnd
	movl	%eax, -188(%rbp)
## BB#73:
	cmpl	$0, -188(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB1_75
## BB#74:
	movl	$1, %eax
	movl	$315, %esi              ## imm = 0x13B
	leaq	L___func__.calc_residual(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-188(%rbp), %r8d
	movq	%rdi, -976(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-976(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB1_97
LBB1_75:
	jmp	LBB1_76
LBB1_76:
	movl	$0, -180(%rbp)
LBB1_77:                                ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB1_79 Depth 2
                                        ##       Child Loop BB1_81 Depth 3
	cmpl	$5, -180(%rbp)
	jge	LBB1_88
## BB#78:                               ##   in Loop: Header=BB1_77 Depth=1
	movl	$0, -184(%rbp)
LBB1_79:                                ##   Parent Loop BB1_77 Depth=1
                                        ## =>  This Loop Header: Depth=2
                                        ##       Child Loop BB1_81 Depth 3
	cmpl	$5, -184(%rbp)
	jge	LBB1_86
## BB#80:                               ##   in Loop: Header=BB1_79 Depth=2
	movl	$0, -176(%rbp)
LBB1_81:                                ##   Parent Loop BB1_77 Depth=1
                                        ##     Parent Loop BB1_79 Depth=2
                                        ## =>    This Inner Loop Header: Depth=3
	cmpl	$2, -176(%rbp)
	jge	LBB1_84
## BB#82:                               ##   in Loop: Header=BB1_81 Depth=3
	movl	$3, %edx
	movl	$2, %ecx
	movq	-16(%rbp), %rdi
	movl	-180(%rbp), %eax
	movl	-184(%rbp), %esi
	movq	%rdi, -984(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$2, %edx
	leaq	_cm(%rip), %r8
	movslq	-176(%rbp), %r9
	movsd	(%r8,%r9,8), %xmm0      ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	%eax, -988(%rbp)        ## 4-byte Spill
	movsd	%xmm0, -1000(%rbp)      ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movsd	%xmm0, -1008(%rbp)      ## 8-byte Spill
	callq	_phi_index
	movl	$2, %edx
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-1008(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm0
	movsd	-1000(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	mulsd	%xmm0, %xmm1
	movq	-984(%rbp), %rdi        ## 8-byte Reload
	movl	-988(%rbp), %esi        ## 4-byte Reload
	movaps	%xmm1, %xmm0
	callq	_VecSetValue
	movl	$3, %edx
	movl	%eax, -188(%rbp)
	movq	-16(%rbp), %rdi
	movl	-180(%rbp), %eax
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %ecx
	movq	%rdi, -1016(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	leaq	_cm(%rip), %r8
	movslq	-176(%rbp), %r9
	movsd	(%r8,%r9,8), %xmm0      ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movl	-176(%rbp), %edx
	movl	%eax, -1020(%rbp)       ## 4-byte Spill
	movsd	%xmm0, -1032(%rbp)      ## 8-byte Spill
	callq	_phi_index
	movl	$2, %edx
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-180(%rbp), %edi
	movl	-184(%rbp), %esi
	movsd	%xmm0, -1040(%rbp)      ## 8-byte Spill
	callq	_phi_index
	movl	$2, %edx
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-1040(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm0
	movsd	-1032(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	mulsd	%xmm0, %xmm1
	movq	-1016(%rbp), %rdi       ## 8-byte Reload
	movl	-1020(%rbp), %esi       ## 4-byte Reload
	movaps	%xmm1, %xmm0
	callq	_VecSetValue
	movl	%eax, -188(%rbp)
## BB#83:                               ##   in Loop: Header=BB1_81 Depth=3
	movl	-176(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -176(%rbp)
	jmp	LBB1_81
LBB1_84:                                ##   in Loop: Header=BB1_79 Depth=2
	jmp	LBB1_85
LBB1_85:                                ##   in Loop: Header=BB1_79 Depth=2
	movl	-184(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -184(%rbp)
	jmp	LBB1_79
LBB1_86:                                ##   in Loop: Header=BB1_77 Depth=1
	jmp	LBB1_87
LBB1_87:                                ##   in Loop: Header=BB1_77 Depth=1
	movl	-180(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -180(%rbp)
	jmp	LBB1_77
LBB1_88:
	movq	-16(%rbp), %rdi
	callq	_VecAssemblyBegin
	movl	%eax, -188(%rbp)
## BB#89:
	cmpl	$0, -188(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB1_91
## BB#90:
	movl	$1, %eax
	movl	$331, %esi              ## imm = 0x14B
	leaq	L___func__.calc_residual(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-188(%rbp), %r8d
	movq	%rdi, -1048(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1048(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB1_97
LBB1_91:
	jmp	LBB1_92
LBB1_92:
	movq	-16(%rbp), %rdi
	callq	_VecAssemblyEnd
	movl	%eax, -188(%rbp)
## BB#93:
	cmpl	$0, -188(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB1_95
## BB#94:
	movl	$1, %eax
	movl	$332, %esi              ## imm = 0x14C
	leaq	L___func__.calc_residual(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-188(%rbp), %r8d
	movq	%rdi, -1056(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1056(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB1_97
LBB1_95:
	jmp	LBB1_96
LBB1_96:
	movl	-188(%rbp), %eax
	movl	%eax, -4(%rbp)
LBB1_97:
	movl	-4(%rbp), %eax
	addq	$1072, %rsp             ## imm = 0x430
	popq	%rbp
	retq
	.cfi_endproc

	.section	__TEXT,__literal8,8byte_literals
	.align	3
LCPI2_0:
	.quad	4607182418800017408     ## double 1
LCPI2_1:
	.quad	4611686018427387904     ## double 2
LCPI2_2:
	.quad	4576918229304087675     ## double 0.01
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_calc_jacobian
	.align	4, 0x90
_calc_jacobian:                         ## @calc_jacobian
	.cfi_startproc
## BB#0:
	pushq	%rbp
Ltmp7:
	.cfi_def_cfa_offset 16
Ltmp8:
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
Ltmp9:
	.cfi_def_cfa_register %rbp
	subq	$1616, %rsp             ## imm = 0x650
	movq	16(%rbp), %rax
	movq	%rdi, -16(%rbp)
	movq	%rsi, -24(%rbp)
	movq	%rdx, -32(%rbp)
	movsd	%xmm0, -40(%rbp)
	movq	%rcx, -48(%rbp)
	movq	%r8, -56(%rbp)
	movq	%r9, -64(%rbp)
	movq	%rax, -72(%rbp)
	movq	-32(%rbp), %rax
	movq	%rax, -80(%rbp)
	movq	-32(%rbp), %rax
	addq	$2400, %rax             ## imm = 0x960
	movq	%rax, -88(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, -96(%rbp)
	movl	$0, -100(%rbp)
	movl	$0, -104(%rbp)
LBB2_1:                                 ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB2_3 Depth 2
                                        ##       Child Loop BB2_5 Depth 3
                                        ##         Child Loop BB2_7 Depth 4
                                        ##         Child Loop BB2_115 Depth 4
	cmpl	$5, -104(%rbp)
	jge	LBB2_132
## BB#2:                                ##   in Loop: Header=BB2_1 Depth=1
	movl	$0, -108(%rbp)
LBB2_3:                                 ##   Parent Loop BB2_1 Depth=1
                                        ## =>  This Loop Header: Depth=2
                                        ##       Child Loop BB2_5 Depth 3
                                        ##         Child Loop BB2_7 Depth 4
                                        ##         Child Loop BB2_115 Depth 4
	cmpl	$5, -108(%rbp)
	jge	LBB2_130
## BB#4:                                ##   in Loop: Header=BB2_3 Depth=2
	movl	$0, -112(%rbp)
LBB2_5:                                 ##   Parent Loop BB2_1 Depth=1
                                        ##     Parent Loop BB2_3 Depth=2
                                        ## =>    This Loop Header: Depth=3
                                        ##         Child Loop BB2_7 Depth 4
                                        ##         Child Loop BB2_115 Depth 4
	cmpl	$3, -112(%rbp)
	jge	LBB2_128
## BB#6:                                ##   in Loop: Header=BB2_5 Depth=3
	movl	$0, -116(%rbp)
LBB2_7:                                 ##   Parent Loop BB2_1 Depth=1
                                        ##     Parent Loop BB2_3 Depth=2
                                        ##       Parent Loop BB2_5 Depth=3
                                        ## =>      This Inner Loop Header: Depth=4
	cmpl	$2, -116(%rbp)
	jge	LBB2_74
## BB#8:                                ##   in Loop: Header=BB2_7 Depth=4
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -128(%rbp)
	movsd	%xmm0, -136(%rbp)
	movsd	%xmm0, -144(%rbp)
	movsd	%xmm0, -152(%rbp)
	movsd	%xmm0, -160(%rbp)
	movsd	%xmm0, -192(%rbp)
	movsd	%xmm0, -200(%rbp)
	cmpl	$4, -104(%rbp)
	jge	LBB2_18
## BB#9:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -248(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -256(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI2_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	-256(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm2
	movsd	-248(%rbp), %xmm3       ## 8-byte Reload
                                        ## xmm3 = mem[0],zero
	mulsd	%xmm2, %xmm3
	divsd	%xmm1, %xmm3
	divsd	%xmm0, %xmm3
	mulsd	-40(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movsd	%xmm3, -128(%rbp)
	movsd	-128(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -264(%rbp)       ## 8-byte Spill
	callq	_c_index
	leaq	_z(%rip), %r8
	movslq	%eax, %r9
	movq	-80(%rbp), %r10
	movsd	-264(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	divsd	(%r10,%r9,8), %xmm0
	movsd	%xmm0, -136(%rbp)
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	mulsd	-128(%rbp), %xmm0
	movsd	%xmm0, -152(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -272(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -276(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-136(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-272(%rbp), %rdi        ## 8-byte Reload
	movl	-276(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#10:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_12
## BB#11:
	movl	$1, %eax
	movl	$375, %esi              ## imm = 0x177
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -288(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-288(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_12:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_13
LBB2_13:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -296(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -300(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-152(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-296(%rbp), %rdi        ## 8-byte Reload
	movl	-300(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#14:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_16
## BB#15:
	movl	$1, %eax
	movl	$378, %esi              ## imm = 0x17A
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -312(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-312(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_16:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_17
LBB2_17:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
LBB2_18:                                ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -104(%rbp)
	jle	LBB2_28
## BB#19:                               ##   in Loop: Header=BB2_7 Depth=4
	movl	-104(%rbp), %eax
	subl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %edi
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %eax
	subl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -320(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -328(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI2_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	-328(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm2
	movsd	-320(%rbp), %xmm3       ## 8-byte Reload
                                        ## xmm3 = mem[0],zero
	mulsd	%xmm2, %xmm3
	divsd	%xmm1, %xmm3
	divsd	%xmm0, %xmm3
	mulsd	-40(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movsd	%xmm3, -128(%rbp)
	movsd	-128(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -336(%rbp)       ## 8-byte Spill
	callq	_c_index
	leaq	_z(%rip), %r8
	movslq	%eax, %r9
	movq	-80(%rbp), %r10
	movsd	-336(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	divsd	(%r10,%r9,8), %xmm0
	movsd	%xmm0, -144(%rbp)
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	mulsd	-128(%rbp), %xmm0
	movsd	%xmm0, -160(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	subl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -344(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -348(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-144(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-344(%rbp), %rdi        ## 8-byte Reload
	movl	-348(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#20:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_22
## BB#21:
	movl	$1, %eax
	movl	$387, %esi              ## imm = 0x183
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -360(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-360(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_22:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_23
LBB2_23:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	subl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -368(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -372(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-160(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-368(%rbp), %rdi        ## 8-byte Reload
	movl	-372(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#24:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_26
## BB#25:
	movl	$1, %eax
	movl	$390, %esi              ## imm = 0x186
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -384(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-384(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_26:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_27
LBB2_27:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
LBB2_28:                                ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$4, -108(%rbp)
	jge	LBB2_38
## BB#29:                               ##   in Loop: Header=BB2_7 Depth=4
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -392(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %eax
	addl	$1, %eax
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -400(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI2_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	-400(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm2
	movsd	-392(%rbp), %xmm3       ## 8-byte Reload
                                        ## xmm3 = mem[0],zero
	mulsd	%xmm2, %xmm3
	divsd	%xmm1, %xmm3
	divsd	%xmm0, %xmm3
	mulsd	-40(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movsd	%xmm3, -168(%rbp)
	movsd	-168(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -408(%rbp)       ## 8-byte Spill
	callq	_c_index
	leaq	_z(%rip), %r8
	movslq	%eax, %r9
	movq	-80(%rbp), %r10
	movsd	-408(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	divsd	(%r10,%r9,8), %xmm0
	movsd	%xmm0, -176(%rbp)
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	mulsd	-168(%rbp), %xmm0
	movsd	%xmm0, -192(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %ecx
	addl	$1, %ecx
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %esi
	movq	%rdi, -416(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%esi, -420(%rbp)        ## 4-byte Spill
	movl	%ecx, %esi
	movl	-420(%rbp), %ecx        ## 4-byte Reload
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -424(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-176(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-416(%rbp), %rdi        ## 8-byte Reload
	movl	-424(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#30:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_32
## BB#31:
	movl	$1, %eax
	movl	$399, %esi              ## imm = 0x18F
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -432(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-432(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_32:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_33
LBB2_33:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %ecx
	addl	$1, %ecx
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %esi
	movq	%rdi, -440(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%esi, -444(%rbp)        ## 4-byte Spill
	movl	%ecx, %esi
	movl	-444(%rbp), %ecx        ## 4-byte Reload
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -448(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-192(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-440(%rbp), %rdi        ## 8-byte Reload
	movl	-448(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#34:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_36
## BB#35:
	movl	$1, %eax
	movl	$402, %esi              ## imm = 0x192
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -456(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-456(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_36:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_37
LBB2_37:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
LBB2_38:                                ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -108(%rbp)
	jle	LBB2_48
## BB#39:                               ##   in Loop: Header=BB2_7 Depth=4
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %eax
	subl	$1, %eax
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %esi
	callq	_c_index
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %eax
	subl	$1, %eax
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -464(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -472(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI2_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	-472(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm2
	movsd	-464(%rbp), %xmm3       ## 8-byte Reload
                                        ## xmm3 = mem[0],zero
	mulsd	%xmm2, %xmm3
	divsd	%xmm1, %xmm3
	divsd	%xmm0, %xmm3
	mulsd	-40(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movsd	%xmm3, -168(%rbp)
	movsd	-168(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -480(%rbp)       ## 8-byte Spill
	callq	_c_index
	leaq	_z(%rip), %r8
	movslq	%eax, %r9
	movq	-80(%rbp), %r10
	movsd	-480(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	divsd	(%r10,%r9,8), %xmm0
	movsd	%xmm0, -184(%rbp)
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	mulsd	-168(%rbp), %xmm0
	movsd	%xmm0, -200(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %ecx
	subl	$1, %ecx
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %esi
	movq	%rdi, -488(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%esi, -492(%rbp)        ## 4-byte Spill
	movl	%ecx, %esi
	movl	-492(%rbp), %ecx        ## 4-byte Reload
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -496(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-184(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-488(%rbp), %rdi        ## 8-byte Reload
	movl	-496(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#40:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_42
## BB#41:
	movl	$1, %eax
	movl	$411, %esi              ## imm = 0x19B
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -504(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-504(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_42:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_43
LBB2_43:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %ecx
	subl	$1, %ecx
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %esi
	movq	%rdi, -512(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%esi, -516(%rbp)        ## 4-byte Spill
	movl	%ecx, %esi
	movl	-516(%rbp), %ecx        ## 4-byte Reload
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -520(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-200(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-512(%rbp), %rdi        ## 8-byte Reload
	movl	-520(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#44:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_46
## BB#45:
	movl	$1, %eax
	movl	$414, %esi              ## imm = 0x19E
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -528(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-528(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_46:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_47
LBB2_47:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
LBB2_48:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	callq	_al_index
	movslq	%eax, %rcx
	movq	-88(%rbp), %r8
	movsd	(%r8,%rcx,8), %xmm0     ## xmm0 = mem[0],zero
	addsd	-136(%rbp), %xmm0
	addsd	-144(%rbp), %xmm0
	addsd	-176(%rbp), %xmm0
	addsd	-184(%rbp), %xmm0
	movsd	%xmm0, -208(%rbp)
	movsd	-152(%rbp), %xmm0       ## xmm0 = mem[0],zero
	addsd	-160(%rbp), %xmm0
	addsd	-192(%rbp), %xmm0
	addsd	-200(%rbp), %xmm0
	movsd	%xmm0, -216(%rbp)
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	1800(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	addsd	-208(%rbp), %xmm0
	movsd	%xmm0, -208(%rbp)
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	movl	$2, %ecx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	5400(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	addsd	-216(%rbp), %xmm0
	movsd	%xmm0, -216(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movq	%rdi, -536(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -540(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, -544(%rbp)        ## 4-byte Spill
	callq	_c_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	1800(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	mulsd	-40(%rbp), %xmm0
	movq	-536(%rbp), %rdi        ## 8-byte Reload
	movl	-540(%rbp), %esi        ## 4-byte Reload
	movl	-544(%rbp), %edx        ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#49:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_51
## BB#50:
	movl	$1, %eax
	movl	$426, %esi              ## imm = 0x1AA
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -552(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-552(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_51:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_52
LBB2_52:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -560(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$2, %ecx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	%eax, -564(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, -568(%rbp)        ## 4-byte Spill
	callq	_c_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	3600(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	movq	-560(%rbp), %rdi        ## 8-byte Reload
	movl	-564(%rbp), %esi        ## 4-byte Reload
	movl	-568(%rbp), %edx        ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#53:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_55
## BB#54:
	movl	$1, %eax
	movl	$429, %esi              ## imm = 0x1AD
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -576(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-576(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_55:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_56
LBB2_56:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	$2, %ecx
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movq	%rdi, -584(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -588(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, -592(%rbp)        ## 4-byte Spill
	callq	_c_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	5400(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	mulsd	-40(%rbp), %xmm0
	movq	-584(%rbp), %rdi        ## 8-byte Reload
	movl	-588(%rbp), %esi        ## 4-byte Reload
	movl	-592(%rbp), %edx        ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#57:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_59
## BB#58:
	movl	$1, %eax
	movl	$432, %esi              ## imm = 0x1B0
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -600(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-600(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_59:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_60
LBB2_60:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -608(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$2, %ecx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	%eax, -612(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, -616(%rbp)        ## 4-byte Spill
	callq	_c_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	5400(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	mulsd	-40(%rbp), %xmm0
	movq	-608(%rbp), %rdi        ## 8-byte Reload
	movl	-612(%rbp), %esi        ## 4-byte Reload
	movl	-616(%rbp), %edx        ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#61:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_63
## BB#62:
	movl	$1, %eax
	movl	$435, %esi              ## imm = 0x1B3
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -624(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-624(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_63:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_64
LBB2_64:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	$2, %ecx
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movq	%rdi, -632(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$4, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -636(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$2, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %ecx
	movl	%eax, -640(%rbp)        ## 4-byte Spill
	callq	_c_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	mulsd	-40(%rbp), %xmm0
	movq	-632(%rbp), %rdi        ## 8-byte Reload
	movl	-636(%rbp), %esi        ## 4-byte Reload
	movl	-640(%rbp), %edx        ## 4-byte Reload
	callq	_MatSetValue
	movl	$2, %ecx
	movl	%eax, -120(%rbp)
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movq	%rdi, -648(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$4, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -652(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, -656(%rbp)        ## 4-byte Spill
	callq	_c_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	movq	-648(%rbp), %rdi        ## 8-byte Reload
	movl	-652(%rbp), %esi        ## 4-byte Reload
	movl	-656(%rbp), %edx        ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -664(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -668(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-208(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-664(%rbp), %rdi        ## 8-byte Reload
	movl	-668(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#65:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_67
## BB#66:
	movl	$1, %eax
	movl	$449, %esi              ## imm = 0x1C1
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -680(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-680(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_67:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_68
LBB2_68:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -688(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -692(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-216(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-688(%rbp), %rdi        ## 8-byte Reload
	movl	-692(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#69:                               ##   in Loop: Header=BB2_7 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_71
## BB#70:
	movl	$1, %eax
	movl	$452, %esi              ## imm = 0x1C4
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -704(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-704(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_71:                                ##   in Loop: Header=BB2_7 Depth=4
	jmp	LBB2_72
LBB2_72:                                ##   in Loop: Header=BB2_7 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
## BB#73:                               ##   in Loop: Header=BB2_7 Depth=4
	movl	-116(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -116(%rbp)
	jmp	LBB2_7
LBB2_74:                                ##   in Loop: Header=BB2_5 Depth=3
	xorps	%xmm0, %xmm0
	movl	$2, -116(%rbp)
	movsd	%xmm0, -128(%rbp)
	movsd	%xmm0, -136(%rbp)
	movsd	%xmm0, -144(%rbp)
	movsd	%xmm0, -152(%rbp)
	movsd	%xmm0, -160(%rbp)
	movsd	%xmm0, -192(%rbp)
	movsd	%xmm0, -200(%rbp)
	cmpl	$4, -104(%rbp)
	jge	LBB2_84
## BB#75:                               ##   in Loop: Header=BB2_5 Depth=3
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -712(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -720(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI2_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	-720(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm2
	movsd	-712(%rbp), %xmm3       ## 8-byte Reload
                                        ## xmm3 = mem[0],zero
	mulsd	%xmm2, %xmm3
	divsd	%xmm1, %xmm3
	divsd	%xmm0, %xmm3
	mulsd	-40(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movsd	%xmm3, -128(%rbp)
	movsd	-128(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -728(%rbp)       ## 8-byte Spill
	callq	_c_index
	leaq	_z(%rip), %r8
	movslq	%eax, %r9
	movq	-80(%rbp), %r10
	movsd	-728(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	divsd	(%r10,%r9,8), %xmm0
	movsd	%xmm0, -136(%rbp)
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	mulsd	-128(%rbp), %xmm0
	movsd	%xmm0, -152(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -736(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -740(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-136(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-736(%rbp), %rdi        ## 8-byte Reload
	movl	-740(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#76:                               ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_78
## BB#77:
	movl	$1, %eax
	movl	$472, %esi              ## imm = 0x1D8
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -752(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-752(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_78:                                ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_79
LBB2_79:                                ##   in Loop: Header=BB2_5 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -760(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -764(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-152(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-760(%rbp), %rdi        ## 8-byte Reload
	movl	-764(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#80:                               ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_82
## BB#81:
	movl	$1, %eax
	movl	$475, %esi              ## imm = 0x1DB
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -776(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-776(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_82:                                ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_83
LBB2_83:                                ##   in Loop: Header=BB2_5 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
LBB2_84:                                ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -104(%rbp)
	jle	LBB2_94
## BB#85:                               ##   in Loop: Header=BB2_5 Depth=3
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %eax
	subl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %edi
	movsd	%xmm0, -784(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -792(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI2_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	-792(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm2
	movsd	-784(%rbp), %xmm3       ## 8-byte Reload
                                        ## xmm3 = mem[0],zero
	mulsd	%xmm2, %xmm3
	divsd	%xmm1, %xmm3
	divsd	%xmm0, %xmm3
	mulsd	-40(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movsd	%xmm3, -128(%rbp)
	movsd	-128(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -800(%rbp)       ## 8-byte Spill
	callq	_c_index
	leaq	_z(%rip), %r8
	movslq	%eax, %r9
	movq	-80(%rbp), %r10
	movsd	-800(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	divsd	(%r10,%r9,8), %xmm0
	movsd	%xmm0, -144(%rbp)
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	mulsd	-128(%rbp), %xmm0
	movsd	%xmm0, -160(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	subl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -808(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -812(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-144(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-808(%rbp), %rdi        ## 8-byte Reload
	movl	-812(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#86:                               ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_88
## BB#87:
	movl	$1, %eax
	movl	$484, %esi              ## imm = 0x1E4
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -824(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-824(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_88:                                ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_89
LBB2_89:                                ##   in Loop: Header=BB2_5 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	subl	$1, %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movq	%rdi, -832(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -836(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-160(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-832(%rbp), %rdi        ## 8-byte Reload
	movl	-836(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#90:                               ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_92
## BB#91:
	movl	$1, %eax
	movl	$487, %esi              ## imm = 0x1E7
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -848(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-848(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_92:                                ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_93
LBB2_93:                                ##   in Loop: Header=BB2_5 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
LBB2_94:                                ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$4, -108(%rbp)
	jge	LBB2_104
## BB#95:                               ##   in Loop: Header=BB2_5 Depth=3
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -856(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %eax
	addl	$1, %eax
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -864(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI2_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	-864(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm2
	movsd	-856(%rbp), %xmm3       ## 8-byte Reload
                                        ## xmm3 = mem[0],zero
	mulsd	%xmm2, %xmm3
	divsd	%xmm1, %xmm3
	divsd	%xmm0, %xmm3
	mulsd	-40(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movsd	%xmm3, -168(%rbp)
	movsd	-168(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -872(%rbp)       ## 8-byte Spill
	callq	_c_index
	leaq	_z(%rip), %r8
	movslq	%eax, %r9
	movq	-80(%rbp), %r10
	movsd	-872(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	divsd	(%r10,%r9,8), %xmm0
	movsd	%xmm0, -176(%rbp)
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	mulsd	-168(%rbp), %xmm0
	movsd	%xmm0, -192(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %ecx
	addl	$1, %ecx
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %esi
	movq	%rdi, -880(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%esi, -884(%rbp)        ## 4-byte Spill
	movl	%ecx, %esi
	movl	-884(%rbp), %ecx        ## 4-byte Reload
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -888(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-176(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-880(%rbp), %rdi        ## 8-byte Reload
	movl	-888(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#96:                               ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_98
## BB#97:
	movl	$1, %eax
	movl	$496, %esi              ## imm = 0x1F0
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -896(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-896(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_98:                                ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_99
LBB2_99:                                ##   in Loop: Header=BB2_5 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %ecx
	addl	$1, %ecx
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %esi
	movq	%rdi, -904(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%esi, -908(%rbp)        ## 4-byte Spill
	movl	%ecx, %esi
	movl	-908(%rbp), %ecx        ## 4-byte Reload
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -912(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-192(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-904(%rbp), %rdi        ## 8-byte Reload
	movl	-912(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#100:                              ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_102
## BB#101:
	movl	$1, %eax
	movl	$499, %esi              ## imm = 0x1F3
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -920(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-920(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_102:                               ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_103
LBB2_103:                               ##   in Loop: Header=BB2_5 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
LBB2_104:                               ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -108(%rbp)
	jle	LBB2_114
## BB#105:                              ##   in Loop: Header=BB2_5 Depth=3
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-48(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %eax
	subl	$1, %eax
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movl	%eax, %esi
	movsd	%xmm0, -928(%rbp)       ## 8-byte Spill
	callq	_c_index
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -936(%rbp)       ## 8-byte Spill
	callq	_c_index
	movsd	LCPI2_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movslq	%eax, %r8
	movq	-96(%rbp), %r9
	movsd	-936(%rbp), %xmm2       ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	addsd	(%r9,%r8,8), %xmm2
	movsd	-928(%rbp), %xmm3       ## 8-byte Reload
                                        ## xmm3 = mem[0],zero
	mulsd	%xmm2, %xmm3
	divsd	%xmm1, %xmm3
	divsd	%xmm0, %xmm3
	mulsd	-40(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movsd	%xmm3, -168(%rbp)
	movsd	-168(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -944(%rbp)       ## 8-byte Spill
	callq	_c_index
	leaq	_z(%rip), %r8
	movslq	%eax, %r9
	movq	-80(%rbp), %r10
	movsd	-944(%rbp), %xmm0       ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	divsd	(%r10,%r9,8), %xmm0
	movsd	%xmm0, -184(%rbp)
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	mulsd	-168(%rbp), %xmm0
	movsd	%xmm0, -200(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %ecx
	subl	$1, %ecx
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %esi
	movq	%rdi, -952(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%esi, -956(%rbp)        ## 4-byte Spill
	movl	%ecx, %esi
	movl	-956(%rbp), %ecx        ## 4-byte Reload
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -960(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-184(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-952(%rbp), %rdi        ## 8-byte Reload
	movl	-960(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#106:                              ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_108
## BB#107:
	movl	$1, %eax
	movl	$508, %esi              ## imm = 0x1FC
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -968(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-968(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_108:                               ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_109
LBB2_109:                               ##   in Loop: Header=BB2_5 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %ecx
	subl	$1, %ecx
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %esi
	movq	%rdi, -976(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%esi, -980(%rbp)        ## 4-byte Spill
	movl	%ecx, %esi
	movl	-980(%rbp), %ecx        ## 4-byte Reload
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -984(%rbp)        ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-200(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-976(%rbp), %rdi        ## 8-byte Reload
	movl	-984(%rbp), %esi        ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#110:                              ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_112
## BB#111:
	movl	$1, %eax
	movl	$511, %esi              ## imm = 0x1FF
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -992(%rbp)        ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-992(%rbp), %r10        ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_112:                               ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_113
LBB2_113:                               ##   in Loop: Header=BB2_5 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
LBB2_114:                               ##   in Loop: Header=BB2_5 Depth=3
	xorl	%edx, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	callq	_al_index
	movl	$1, %edx
	movsd	LCPI2_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %rcx
	movq	-88(%rbp), %r8
	subsd	(%r8,%rcx,8), %xmm0
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movsd	%xmm0, -1000(%rbp)      ## 8-byte Spill
	callq	_al_index
	movslq	%eax, %rcx
	movq	-88(%rbp), %r8
	movsd	-1000(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r8,%rcx,8), %xmm0
	addsd	-136(%rbp), %xmm0
	addsd	-144(%rbp), %xmm0
	addsd	-176(%rbp), %xmm0
	addsd	-184(%rbp), %xmm0
	movsd	%xmm0, -208(%rbp)
	movsd	-152(%rbp), %xmm0       ## xmm0 = mem[0],zero
	addsd	-160(%rbp), %xmm0
	addsd	-192(%rbp), %xmm0
	addsd	-200(%rbp), %xmm0
	movsd	%xmm0, -216(%rbp)
	movl	$0, -116(%rbp)
LBB2_115:                               ##   Parent Loop BB2_1 Depth=1
                                        ##     Parent Loop BB2_3 Depth=2
                                        ##       Parent Loop BB2_5 Depth=3
                                        ## =>      This Inner Loop Header: Depth=4
	cmpl	$2, -116(%rbp)
	jge	LBB2_118
## BB#116:                              ##   in Loop: Header=BB2_115 Depth=4
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	3600(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	movsd	-208(%rbp), %xmm1       ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movsd	%xmm1, -208(%rbp)
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	-112(%rbp), %ecx
	callq	_c_index
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	5400(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	addsd	-216(%rbp), %xmm0
	movsd	%xmm0, -216(%rbp)
## BB#117:                              ##   in Loop: Header=BB2_115 Depth=4
	movl	-116(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -116(%rbp)
	jmp	LBB2_115
LBB2_118:                               ##   in Loop: Header=BB2_5 Depth=3
	movl	$2, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %ecx
	callq	_c_index
	movl	$2, %edx
	shll	$1, %eax
	movslq	%eax, %r8
	movq	-56(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movl	%edx, -1004(%rbp)       ## 4-byte Spill
	callq	_pow
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %ecx
	movl	-1004(%rbp), %edx       ## 4-byte Reload
	movsd	%xmm0, -1016(%rbp)      ## 8-byte Spill
	callq	_c_index
	movl	$2, %edx
	shll	$1, %eax
	addl	$1, %eax
	movslq	%eax, %r8
	movq	-56(%rbp), %r9
	movsd	(%r9,%r8,8), %xmm0      ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movl	%edx, -1020(%rbp)       ## 4-byte Spill
	callq	_pow
	movsd	-1016(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	%xmm0, %xmm1
	sqrtsd	%xmm1, %xmm0
	movsd	%xmm0, -128(%rbp)
	movsd	-128(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %ecx
	movl	-1020(%rbp), %edx       ## 4-byte Reload
	movsd	%xmm0, -1032(%rbp)      ## 8-byte Spill
	callq	_c_index
	movl	$2, %edx
	leaq	_cbath(%rip), %r8
	movslq	%eax, %r9
	movq	-96(%rbp), %r10
	movsd	(%r10,%r9,8), %xmm0     ## xmm0 = mem[0],zero
	movslq	-112(%rbp), %r9
	addsd	(%r8,%r9,8), %xmm0
	movsd	-1032(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	mulsd	%xmm0, %xmm1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %ecx
	movsd	%xmm1, -1040(%rbp)      ## 8-byte Spill
	callq	_c_index
	movl	$2, %edx
	movsd	LCPI2_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-80(%rbp), %r9
	mulsd	(%r9,%r8,8), %xmm0
	movsd	-1040(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	divsd	%xmm0, %xmm1
	mulsd	-40(%rbp), %xmm1
	movsd	-208(%rbp), %xmm0       ## xmm0 = mem[0],zero
	subsd	%xmm1, %xmm0
	movsd	%xmm0, -208(%rbp)
	movsd	-128(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %ecx
	movsd	%xmm0, -1048(%rbp)      ## 8-byte Spill
	callq	_c_index
	movl	$2, %ecx
	movsd	LCPI2_1(%rip), %xmm0    ## xmm0 = mem[0],zero
	leaq	_z(%rip), %r8
	leaq	_cbath(%rip), %r9
	movslq	%eax, %r10
	movq	-96(%rbp), %r11
	movsd	(%r11,%r10,8), %xmm1    ## xmm1 = mem[0],zero
	movslq	-112(%rbp), %r10
	addsd	(%r9,%r10,8), %xmm1
	movsd	-1048(%rbp), %xmm2      ## 8-byte Reload
                                        ## xmm2 = mem[0],zero
	mulsd	%xmm1, %xmm2
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm1
	mulsd	%xmm1, %xmm2
	divsd	%xmm0, %xmm2
	mulsd	-40(%rbp), %xmm2
	movsd	-216(%rbp), %xmm0       ## xmm0 = mem[0],zero
	subsd	%xmm2, %xmm0
	movsd	%xmm0, -216(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movq	%rdi, -1056(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$2, %ecx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	%eax, -1060(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-208(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-1056(%rbp), %rdi       ## 8-byte Reload
	movl	-1060(%rbp), %esi       ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#119:                              ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_121
## BB#120:
	movl	$1, %eax
	movl	$532, %esi              ## imm = 0x214
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1072(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1072(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_121:                               ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_122
LBB2_122:                               ##   in Loop: Header=BB2_5 Depth=3
	movl	$2, %ecx
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movq	%rdi, -1080(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	$2, %ecx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	%eax, -1084(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-216(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-1080(%rbp), %rdi       ## 8-byte Reload
	movl	-1084(%rbp), %esi       ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#123:                              ##   in Loop: Header=BB2_5 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_125
## BB#124:
	movl	$1, %eax
	movl	$535, %esi              ## imm = 0x217
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1096(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1096(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_125:                               ##   in Loop: Header=BB2_5 Depth=3
	jmp	LBB2_126
LBB2_126:                               ##   in Loop: Header=BB2_5 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
## BB#127:                              ##   in Loop: Header=BB2_5 Depth=3
	movl	-112(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -112(%rbp)
	jmp	LBB2_5
LBB2_128:                               ##   in Loop: Header=BB2_3 Depth=2
	jmp	LBB2_129
LBB2_129:                               ##   in Loop: Header=BB2_3 Depth=2
	movl	-108(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -108(%rbp)
	jmp	LBB2_3
LBB2_130:                               ##   in Loop: Header=BB2_1 Depth=1
	jmp	LBB2_131
LBB2_131:                               ##   in Loop: Header=BB2_1 Depth=1
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -104(%rbp)
	jmp	LBB2_1
LBB2_132:
	movl	$0, -104(%rbp)
LBB2_133:                               ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB2_135 Depth 2
                                        ##       Child Loop BB2_137 Depth 3
                                        ##         Child Loop BB2_139 Depth 4
                                        ##       Child Loop BB2_153 Depth 3
                                        ##       Child Loop BB2_161 Depth 3
	cmpl	$5, -104(%rbp)
	jge	LBB2_188
## BB#134:                              ##   in Loop: Header=BB2_133 Depth=1
	movl	$0, -108(%rbp)
LBB2_135:                               ##   Parent Loop BB2_133 Depth=1
                                        ## =>  This Loop Header: Depth=2
                                        ##       Child Loop BB2_137 Depth 3
                                        ##         Child Loop BB2_139 Depth 4
                                        ##       Child Loop BB2_153 Depth 3
                                        ##       Child Loop BB2_161 Depth 3
	cmpl	$5, -108(%rbp)
	jge	LBB2_186
## BB#136:                              ##   in Loop: Header=BB2_135 Depth=2
	movl	$0, -112(%rbp)
LBB2_137:                               ##   Parent Loop BB2_133 Depth=1
                                        ##     Parent Loop BB2_135 Depth=2
                                        ## =>    This Loop Header: Depth=3
                                        ##         Child Loop BB2_139 Depth 4
	cmpl	$3, -112(%rbp)
	jge	LBB2_152
## BB#138:                              ##   in Loop: Header=BB2_137 Depth=3
	movl	$0, -116(%rbp)
LBB2_139:                               ##   Parent Loop BB2_133 Depth=1
                                        ##     Parent Loop BB2_135 Depth=2
                                        ##       Parent Loop BB2_137 Depth=3
                                        ## =>      This Inner Loop Header: Depth=4
	cmpl	$2, -116(%rbp)
	jge	LBB2_146
## BB#140:                              ##   in Loop: Header=BB2_139 Depth=4
	movl	$3, %edx
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1104(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -1108(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	leaq	_z(%rip), %r8
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	%eax, -1112(%rbp)       ## 4-byte Spill
	movsd	%xmm0, -1120(%rbp)      ## 8-byte Spill
	callq	_al_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-1120(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	(%r9,%r8,8), %xmm0
	movq	-1104(%rbp), %rdi       ## 8-byte Reload
	movl	-1108(%rbp), %esi       ## 4-byte Reload
	movl	-1112(%rbp), %edx       ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#141:                              ##   in Loop: Header=BB2_139 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_143
## BB#142:
	movl	$1, %eax
	movl	$551, %esi              ## imm = 0x227
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1128(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1128(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_143:                               ##   in Loop: Header=BB2_139 Depth=4
	jmp	LBB2_144
LBB2_144:                               ##   in Loop: Header=BB2_139 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
## BB#145:                              ##   in Loop: Header=BB2_139 Depth=4
	movl	-116(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -116(%rbp)
	jmp	LBB2_139
LBB2_146:                               ##   in Loop: Header=BB2_137 Depth=3
	movl	$3, %edx
	movl	$2, -116(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1136(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-112(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -1140(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	xorl	%edx, %edx
	leaq	_z(%rip), %r8
	movslq	-112(%rbp), %r9
	cvtsi2sdl	(%r8,%r9,4), %xmm0
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	%eax, -1144(%rbp)       ## 4-byte Spill
	movsd	%xmm0, -1152(%rbp)      ## 8-byte Spill
	callq	_al_index
	movl	$1, %edx
	movsd	LCPI2_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	subsd	(%r9,%r8,8), %xmm0
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movsd	%xmm0, -1160(%rbp)      ## 8-byte Spill
	callq	_al_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-1160(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm0
	movsd	-1152(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	mulsd	%xmm0, %xmm1
	movq	-1136(%rbp), %rdi       ## 8-byte Reload
	movl	-1140(%rbp), %esi       ## 4-byte Reload
	movl	-1144(%rbp), %edx       ## 4-byte Reload
	movaps	%xmm1, %xmm0
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#147:                              ##   in Loop: Header=BB2_137 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_149
## BB#148:
	movl	$1, %eax
	movl	$556, %esi              ## imm = 0x22C
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1168(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1168(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_149:                               ##   in Loop: Header=BB2_137 Depth=3
	jmp	LBB2_150
LBB2_150:                               ##   in Loop: Header=BB2_137 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
## BB#151:                              ##   in Loop: Header=BB2_137 Depth=3
	movl	-112(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -112(%rbp)
	jmp	LBB2_137
LBB2_152:                               ##   in Loop: Header=BB2_135 Depth=2
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -216(%rbp)
	movl	$0, -116(%rbp)
LBB2_153:                               ##   Parent Loop BB2_133 Depth=1
                                        ##     Parent Loop BB2_135 Depth=2
                                        ## =>    This Inner Loop Header: Depth=3
	cmpl	$2, -116(%rbp)
	jge	LBB2_156
## BB#154:                              ##   in Loop: Header=BB2_153 Depth=3
	leaq	_cm(%rip), %rax
	movslq	-116(%rbp), %rcx
	movsd	(%rax,%rcx,8), %xmm0    ## xmm0 = mem[0],zero
	movsd	-216(%rbp), %xmm1       ## xmm1 = mem[0],zero
	subsd	%xmm0, %xmm1
	movsd	%xmm1, -216(%rbp)
## BB#155:                              ##   in Loop: Header=BB2_153 Depth=3
	movl	-116(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -116(%rbp)
	jmp	LBB2_153
LBB2_156:                               ##   in Loop: Header=BB2_135 Depth=2
	movl	$3, %edx
	movl	$2, %ecx
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movq	%rdi, -1176(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	$2, %ecx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	%eax, -1180(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-216(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-1176(%rbp), %rdi       ## 8-byte Reload
	movl	-1180(%rbp), %esi       ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#157:                              ##   in Loop: Header=BB2_135 Depth=2
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_159
## BB#158:
	movl	$1, %eax
	movl	$566, %esi              ## imm = 0x236
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1192(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1192(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_159:                               ##   in Loop: Header=BB2_135 Depth=2
	jmp	LBB2_160
LBB2_160:                               ##   in Loop: Header=BB2_135 Depth=2
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movl	$0, -116(%rbp)
LBB2_161:                               ##   Parent Loop BB2_133 Depth=1
                                        ##     Parent Loop BB2_135 Depth=2
                                        ## =>    This Inner Loop Header: Depth=3
	cmpl	$2, -116(%rbp)
	jge	LBB2_184
## BB#162:                              ##   in Loop: Header=BB2_161 Depth=3
	movl	$3, %edx
	movl	$2, %ecx
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movq	%rdi, -1200(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -1204(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	leaq	_cm(%rip), %r8
	movslq	-116(%rbp), %r9
	movsd	(%r8,%r9,8), %xmm0      ## xmm0 = mem[0],zero
	movq	-1200(%rbp), %rdi       ## 8-byte Reload
	movl	-1204(%rbp), %esi       ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#163:                              ##   in Loop: Header=BB2_161 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_165
## BB#164:
	movl	$1, %eax
	movl	$571, %esi              ## imm = 0x23B
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1216(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1216(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_165:                               ##   in Loop: Header=BB2_161 Depth=3
	jmp	LBB2_166
LBB2_166:                               ##   in Loop: Header=BB2_161 Depth=3
	movl	$3, %edx
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1224(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	$2, %ecx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	%eax, -1228(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	leaq	_cm(%rip), %r8
	movslq	-116(%rbp), %r9
	movsd	(%r8,%r9,8), %xmm0      ## xmm0 = mem[0],zero
	movq	-1224(%rbp), %rdi       ## 8-byte Reload
	movl	-1228(%rbp), %esi       ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#167:                              ##   in Loop: Header=BB2_161 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_169
## BB#168:
	movl	$1, %eax
	movl	$574, %esi              ## imm = 0x23E
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1240(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1240(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_169:                               ##   in Loop: Header=BB2_161 Depth=3
	jmp	LBB2_170
LBB2_170:                               ##   in Loop: Header=BB2_161 Depth=3
	movl	$3, %edx
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1248(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$3, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -1252(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	leaq	_cm(%rip), %r8
	movslq	-116(%rbp), %r9
	movsd	(%r8,%r9,8), %xmm0      ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	movq	-1248(%rbp), %rdi       ## 8-byte Reload
	movl	-1252(%rbp), %esi       ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#171:                              ##   in Loop: Header=BB2_161 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_173
## BB#172:
	movl	$1, %eax
	movl	$577, %esi              ## imm = 0x241
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1264(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1264(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_173:                               ##   in Loop: Header=BB2_161 Depth=3
	jmp	LBB2_174
LBB2_174:                               ##   in Loop: Header=BB2_161 Depth=3
	movl	$3, %edx
	movl	$2, %ecx
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movq	%rdi, -1272(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$4, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -1276(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	leaq	_z(%rip), %rsi
	movl	$2, %r8d
	movq	-80(%rbp), %rdi
	movl	-104(%rbp), %edx
	movl	-108(%rbp), %ecx
	movl	%eax, -1280(%rbp)       ## 4-byte Spill
	callq	_cz
	movl	$1, %ecx
	movd	%xmm0, %rsi
	movabsq	$-9223372036854775808, %rdi ## imm = 0x8000000000000000
	xorq	%rdi, %rsi
	movd	%rsi, %xmm0
	movq	-1272(%rbp), %rdi       ## 8-byte Reload
	movl	-1276(%rbp), %esi       ## 4-byte Reload
	movl	-1280(%rbp), %edx       ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#175:                              ##   in Loop: Header=BB2_161 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_177
## BB#176:
	movl	$1, %eax
	movl	$580, %esi              ## imm = 0x244
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1288(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1288(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_177:                               ##   in Loop: Header=BB2_161 Depth=3
	jmp	LBB2_178
LBB2_178:                               ##   in Loop: Header=BB2_161 Depth=3
	movl	$3, %edx
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1296(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$4, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -1300(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	leaq	_z(%rip), %rsi
	movq	-80(%rbp), %rdi
	movl	-104(%rbp), %edx
	movl	-108(%rbp), %ecx
	movl	-116(%rbp), %r8d
	movl	%eax, -1304(%rbp)       ## 4-byte Spill
	callq	_cz
	movl	$1, %ecx
	movq	-1296(%rbp), %rdi       ## 8-byte Reload
	movl	-1300(%rbp), %esi       ## 4-byte Reload
	movl	-1304(%rbp), %edx       ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#179:                              ##   in Loop: Header=BB2_161 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_181
## BB#180:
	movl	$1, %eax
	movl	$583, %esi              ## imm = 0x247
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1312(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1312(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_181:                               ##   in Loop: Header=BB2_161 Depth=3
	jmp	LBB2_182
LBB2_182:                               ##   in Loop: Header=BB2_161 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
## BB#183:                              ##   in Loop: Header=BB2_161 Depth=3
	movl	-116(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -116(%rbp)
	jmp	LBB2_161
LBB2_184:                               ##   in Loop: Header=BB2_135 Depth=2
	jmp	LBB2_185
LBB2_185:                               ##   in Loop: Header=BB2_135 Depth=2
	movl	-108(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -108(%rbp)
	jmp	LBB2_135
LBB2_186:                               ##   in Loop: Header=BB2_133 Depth=1
	jmp	LBB2_187
LBB2_187:                               ##   in Loop: Header=BB2_133 Depth=1
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -104(%rbp)
	jmp	LBB2_133
LBB2_188:
	movl	$0, -104(%rbp)
LBB2_189:                               ## =>This Loop Header: Depth=1
                                        ##     Child Loop BB2_191 Depth 2
                                        ##       Child Loop BB2_193 Depth 3
                                        ##         Child Loop BB2_199 Depth 4
                                        ##         Child Loop BB2_207 Depth 4
                                        ##         Child Loop BB2_215 Depth 4
	cmpl	$5, -104(%rbp)
	jge	LBB2_232
## BB#190:                              ##   in Loop: Header=BB2_189 Depth=1
	movl	$0, -108(%rbp)
LBB2_191:                               ##   Parent Loop BB2_189 Depth=1
                                        ## =>  This Loop Header: Depth=2
                                        ##       Child Loop BB2_193 Depth 3
                                        ##         Child Loop BB2_199 Depth 4
                                        ##         Child Loop BB2_207 Depth 4
                                        ##         Child Loop BB2_215 Depth 4
	cmpl	$5, -108(%rbp)
	jge	LBB2_230
## BB#192:                              ##   in Loop: Header=BB2_191 Depth=2
	movl	$0, -116(%rbp)
LBB2_193:                               ##   Parent Loop BB2_189 Depth=1
                                        ##     Parent Loop BB2_191 Depth=2
                                        ## =>    This Loop Header: Depth=3
                                        ##         Child Loop BB2_199 Depth 4
                                        ##         Child Loop BB2_207 Depth 4
                                        ##         Child Loop BB2_215 Depth 4
	cmpl	$2, -116(%rbp)
	jge	LBB2_228
## BB#194:                              ##   in Loop: Header=BB2_193 Depth=3
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	callq	_al_index
	xorl	%edx, %edx
	movl	$2, %esi
	movslq	%eax, %rcx
	movq	-64(%rbp), %r8
	movsd	7600(%r8,%rcx,8), %xmm0 ## xmm0 = mem[0],zero
	movl	%edx, %edi
	movl	%esi, -1316(%rbp)       ## 4-byte Spill
	movl	%edx, %esi
	movl	-1316(%rbp), %edx       ## 4-byte Reload
	movsd	%xmm0, -1328(%rbp)      ## 8-byte Spill
	callq	_phi_index
	xorl	%edx, %edx
	movslq	%eax, %rcx
	movq	-72(%rbp), %r8
	movsd	40(%r8,%rcx,8), %xmm0   ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movsd	%xmm0, -1336(%rbp)      ## 8-byte Spill
	callq	_al_index
	movl	$1, %edx
	movsd	LCPI2_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %rcx
	movq	-88(%rbp), %r8
	subsd	(%r8,%rcx,8), %xmm0
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movsd	%xmm0, -1344(%rbp)      ## 8-byte Spill
	callq	_al_index
	xorl	%edx, %edx
	movslq	%eax, %rcx
	movq	-88(%rbp), %r8
	movsd	-1344(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r8,%rcx,8), %xmm0
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movl	%edx, -1348(%rbp)       ## 4-byte Spill
	callq	_pow
	movsd	-1336(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	divsd	%xmm0, %xmm1
	movl	-116(%rbp), %edx
	movl	-1348(%rbp), %edi       ## 4-byte Reload
	movl	-1348(%rbp), %esi       ## 4-byte Reload
	movsd	%xmm1, -1360(%rbp)      ## 8-byte Spill
	callq	_phi_index
	movslq	%eax, %rcx
	movq	-72(%rbp), %r8
	movsd	40(%r8,%rcx,8), %xmm0   ## xmm0 = mem[0],zero
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movsd	%xmm0, -1368(%rbp)      ## 8-byte Spill
	callq	_al_index
	movslq	%eax, %rcx
	movq	-88(%rbp), %r8
	movsd	(%r8,%rcx,8), %xmm0     ## xmm0 = mem[0],zero
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	callq	_pow
	movsd	-1368(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	divsd	%xmm0, %xmm1
	movsd	-1360(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	addsd	%xmm1, %xmm0
	movsd	-1328(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	mulsd	%xmm0, %xmm1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movsd	%xmm1, -1376(%rbp)      ## 8-byte Spill
	callq	_al_index
	movl	$4, %edx
	movsd	LCPI2_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %rcx
	movq	-64(%rbp), %r8
	movsd	-1376(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	addsd	8000(%r8,%rcx,8), %xmm1
	mulsd	-40(%rbp), %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -208(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1384(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$4, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movl	%eax, -1388(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	$1, %ecx
	movsd	-208(%rbp), %xmm0       ## xmm0 = mem[0],zero
	movq	-1384(%rbp), %rdi       ## 8-byte Reload
	movl	-1388(%rbp), %esi       ## 4-byte Reload
	movl	%eax, %edx
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#195:                              ##   in Loop: Header=BB2_193 Depth=3
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_197
## BB#196:
	movl	$1, %eax
	movl	$598, %esi              ## imm = 0x256
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1400(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1400(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_197:                               ##   in Loop: Header=BB2_193 Depth=3
	jmp	LBB2_198
LBB2_198:                               ##   in Loop: Header=BB2_193 Depth=3
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movl	$0, -220(%rbp)
LBB2_199:                               ##   Parent Loop BB2_189 Depth=1
                                        ##     Parent Loop BB2_191 Depth=2
                                        ##       Parent Loop BB2_193 Depth=3
                                        ## =>      This Inner Loop Header: Depth=4
	movl	-220(%rbp), %eax
	cmpl	-116(%rbp), %eax
	jge	LBB2_206
## BB#200:                              ##   in Loop: Header=BB2_199 Depth=4
	movl	$4, %edx
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1408(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$4, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-220(%rbp), %ecx
	movl	%eax, -1412(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	%eax, -1416(%rbp)       ## 4-byte Spill
	callq	_al_index
	xorl	%ecx, %ecx
	movl	$2, %edx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	7600(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	movl	%ecx, %edi
	movl	%ecx, %esi
	movsd	%xmm0, -1424(%rbp)      ## 8-byte Spill
	callq	_phi_index
	xorl	%edx, %edx
	movslq	%eax, %r8
	movq	-72(%rbp), %r9
	movsd	-1424(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	40(%r9,%r8,8), %xmm0
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movsd	%xmm0, -1432(%rbp)      ## 8-byte Spill
	callq	_al_index
	movl	$1, %edx
	movsd	LCPI2_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	subsd	(%r9,%r8,8), %xmm0
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movsd	%xmm0, -1440(%rbp)      ## 8-byte Spill
	callq	_al_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-1440(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm0
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movl	%ecx, -1444(%rbp)       ## 4-byte Spill
	callq	_pow
	movsd	-1432(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	divsd	%xmm0, %xmm1
	mulsd	-40(%rbp), %xmm1
	movq	-1408(%rbp), %rdi       ## 8-byte Reload
	movl	-1412(%rbp), %esi       ## 4-byte Reload
	movl	-1416(%rbp), %edx       ## 4-byte Reload
	movaps	%xmm1, %xmm0
	movl	-1444(%rbp), %ecx       ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#201:                              ##   in Loop: Header=BB2_199 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_203
## BB#202:
	movl	$1, %eax
	movl	$603, %esi              ## imm = 0x25B
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1456(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1456(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_203:                               ##   in Loop: Header=BB2_199 Depth=4
	jmp	LBB2_204
LBB2_204:                               ##   in Loop: Header=BB2_199 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
## BB#205:                              ##   in Loop: Header=BB2_199 Depth=4
	movl	-220(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -220(%rbp)
	jmp	LBB2_199
LBB2_206:                               ##   in Loop: Header=BB2_193 Depth=3
	movl	-116(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -224(%rbp)
LBB2_207:                               ##   Parent Loop BB2_189 Depth=1
                                        ##     Parent Loop BB2_191 Depth=2
                                        ##       Parent Loop BB2_193 Depth=3
                                        ## =>      This Inner Loop Header: Depth=4
	cmpl	$2, -224(%rbp)
	jge	LBB2_214
## BB#208:                              ##   in Loop: Header=BB2_207 Depth=4
	movl	$4, %edx
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1464(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$4, %edx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-224(%rbp), %ecx
	movl	%eax, -1468(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	%eax, -1472(%rbp)       ## 4-byte Spill
	callq	_al_index
	xorl	%ecx, %ecx
	movl	$2, %edx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	7600(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	movl	%ecx, %edi
	movl	%ecx, %esi
	movsd	%xmm0, -1480(%rbp)      ## 8-byte Spill
	callq	_phi_index
	xorl	%edx, %edx
	movslq	%eax, %r8
	movq	-72(%rbp), %r9
	movsd	-1480(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	mulsd	40(%r9,%r8,8), %xmm0
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movsd	%xmm0, -1488(%rbp)      ## 8-byte Spill
	callq	_al_index
	movl	$1, %edx
	movsd	LCPI2_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	subsd	(%r9,%r8,8), %xmm0
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movsd	%xmm0, -1496(%rbp)      ## 8-byte Spill
	callq	_al_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-88(%rbp), %r9
	movsd	-1496(%rbp), %xmm0      ## 8-byte Reload
                                        ## xmm0 = mem[0],zero
	subsd	(%r9,%r8,8), %xmm0
	movsd	LCPI2_1(%rip), %xmm1    ## xmm1 = mem[0],zero
	movl	%ecx, -1500(%rbp)       ## 4-byte Spill
	callq	_pow
	movsd	-1488(%rbp), %xmm1      ## 8-byte Reload
                                        ## xmm1 = mem[0],zero
	divsd	%xmm0, %xmm1
	mulsd	-40(%rbp), %xmm1
	movq	-1464(%rbp), %rdi       ## 8-byte Reload
	movl	-1468(%rbp), %esi       ## 4-byte Reload
	movl	-1472(%rbp), %edx       ## 4-byte Reload
	movaps	%xmm1, %xmm0
	movl	-1500(%rbp), %ecx       ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#209:                              ##   in Loop: Header=BB2_207 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_211
## BB#210:
	movl	$1, %eax
	movl	$608, %esi              ## imm = 0x260
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1512(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1512(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_211:                               ##   in Loop: Header=BB2_207 Depth=4
	jmp	LBB2_212
LBB2_212:                               ##   in Loop: Header=BB2_207 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
## BB#213:                              ##   in Loop: Header=BB2_207 Depth=4
	movl	-224(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -224(%rbp)
	jmp	LBB2_207
LBB2_214:                               ##   in Loop: Header=BB2_193 Depth=3
	movl	$0, -228(%rbp)
LBB2_215:                               ##   Parent Loop BB2_189 Depth=1
                                        ##     Parent Loop BB2_191 Depth=2
                                        ##       Parent Loop BB2_193 Depth=3
                                        ## =>      This Inner Loop Header: Depth=4
	cmpl	$3, -228(%rbp)
	jge	LBB2_226
## BB#216:                              ##   in Loop: Header=BB2_215 Depth=4
	movl	$4, %edx
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1520(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	$2, %ecx
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-228(%rbp), %edx
	movl	%eax, -1524(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	%eax, -1528(%rbp)       ## 4-byte Spill
	callq	_al_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	7600(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm0
	movq	-1520(%rbp), %rdi       ## 8-byte Reload
	movl	-1524(%rbp), %esi       ## 4-byte Reload
	movl	-1528(%rbp), %edx       ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#217:                              ##   in Loop: Header=BB2_215 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_219
## BB#218:
	movl	$1, %eax
	movl	$614, %esi              ## imm = 0x266
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1536(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1536(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_219:                               ##   in Loop: Header=BB2_215 Depth=4
	jmp	LBB2_220
LBB2_220:                               ##   in Loop: Header=BB2_215 Depth=4
	movl	$4, %edx
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
	movq	-16(%rbp), %rdi
	movl	-104(%rbp), %eax
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %ecx
	movq	%rdi, -1544(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-228(%rbp), %edx
	movl	-116(%rbp), %ecx
	movl	%eax, -1548(%rbp)       ## 4-byte Spill
	callq	_Ind_1
	movl	-104(%rbp), %edi
	movl	-108(%rbp), %esi
	movl	-116(%rbp), %edx
	movl	%eax, -1552(%rbp)       ## 4-byte Spill
	callq	_al_index
	movl	$1, %ecx
	movslq	%eax, %r8
	movq	-64(%rbp), %r9
	movsd	7600(%r9,%r8,8), %xmm0  ## xmm0 = mem[0],zero
	movd	%xmm0, %r8
	movabsq	$-9223372036854775808, %r9 ## imm = 0x8000000000000000
	xorq	%r9, %r8
	movd	%r8, %xmm0
	mulsd	-40(%rbp), %xmm0
	movq	-1544(%rbp), %rdi       ## 8-byte Reload
	movl	-1548(%rbp), %esi       ## 4-byte Reload
	movl	-1552(%rbp), %edx       ## 4-byte Reload
	callq	_MatSetValue
	movl	%eax, -120(%rbp)
## BB#221:                              ##   in Loop: Header=BB2_215 Depth=4
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_223
## BB#222:
	movl	$1, %eax
	movl	$617, %esi              ## imm = 0x269
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1560(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1560(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_223:                               ##   in Loop: Header=BB2_215 Depth=4
	jmp	LBB2_224
LBB2_224:                               ##   in Loop: Header=BB2_215 Depth=4
	movl	-100(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -100(%rbp)
## BB#225:                              ##   in Loop: Header=BB2_215 Depth=4
	movl	-228(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -228(%rbp)
	jmp	LBB2_215
LBB2_226:                               ##   in Loop: Header=BB2_193 Depth=3
	jmp	LBB2_227
LBB2_227:                               ##   in Loop: Header=BB2_193 Depth=3
	movl	-116(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -116(%rbp)
	jmp	LBB2_193
LBB2_228:                               ##   in Loop: Header=BB2_191 Depth=2
	jmp	LBB2_229
LBB2_229:                               ##   in Loop: Header=BB2_191 Depth=2
	movl	-108(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -108(%rbp)
	jmp	LBB2_191
LBB2_230:                               ##   in Loop: Header=BB2_189 Depth=1
	jmp	LBB2_231
LBB2_231:                               ##   in Loop: Header=BB2_189 Depth=1
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -104(%rbp)
	jmp	LBB2_189
LBB2_232:
	xorl	%esi, %esi
	movq	-16(%rbp), %rdi
	callq	_MatAssemblyBegin
	movl	%eax, -120(%rbp)
## BB#233:
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_235
## BB#234:
	movl	$1, %eax
	movl	$623, %esi              ## imm = 0x26F
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1568(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1568(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_235:
	jmp	LBB2_236
LBB2_236:
	xorl	%esi, %esi
	movq	-16(%rbp), %rdi
	callq	_MatAssemblyEnd
	movl	%eax, -120(%rbp)
## BB#237:
	cmpl	$0, -120(%rbp)
	setne	%al
	xorb	$-1, %al
	xorb	$-1, %al
	andb	$1, %al
	movzbl	%al, %ecx
	movslq	%ecx, %rdx
	cmpq	$0, %rdx
	je	LBB2_239
## BB#238:
	movl	$1, %eax
	movl	$624, %esi              ## imm = 0x270
	leaq	L___func__.calc_jacobian(%rip), %rdx
	leaq	L_.str.20(%rip), %rcx
	leaq	L_.str.21(%rip), %rdi
	movl	-120(%rbp), %r8d
	movq	%rdi, -1576(%rbp)       ## 8-byte Spill
	movl	%eax, %edi
	movl	%eax, %r9d
	movq	-1576(%rbp), %r10       ## 8-byte Reload
	movq	%r10, (%rsp)
	movb	$0, %al
	callq	_PetscError
	movl	%eax, -4(%rbp)
	jmp	LBB2_241
LBB2_239:
	jmp	LBB2_240
LBB2_240:
	leaq	L_.str.24(%rip), %rdi
	movl	$3690, %esi             ## imm = 0xE6A
	movl	-100(%rbp), %edx
	movb	$0, %al
	callq	_printf
	movl	$1, %edi
	movq	-16(%rbp), %rcx
	movl	%eax, -1580(%rbp)       ## 4-byte Spill
	movq	%rcx, -1592(%rbp)       ## 8-byte Spill
	callq	_PETSC_VIEWER_STDOUT_
	movq	-1592(%rbp), %rdi       ## 8-byte Reload
	movq	%rax, %rsi
	callq	_MatView
	leaq	L_.str.25(%rip), %rsi
	leaq	-240(%rbp), %rdx
	movq	_PETSC_COMM_WORLD@GOTPCREL(%rip), %rcx
	movl	(%rcx), %edi
	movl	%eax, -1596(%rbp)       ## 4-byte Spill
	callq	_PetscViewerASCIIOpen
	movq	-16(%rbp), %rdi
	movq	-240(%rbp), %rsi
	movl	%eax, -1600(%rbp)       ## 4-byte Spill
	callq	_MatView
	movl	-120(%rbp), %r8d
	movl	%r8d, -4(%rbp)
	movl	%eax, -1604(%rbp)       ## 4-byte Spill
LBB2_241:
	movl	-4(%rbp), %eax
	addq	$1616, %rsp             ## imm = 0x650
	popq	%rbp
	retq
	.cfi_endproc

	.align	4, 0x90
_VecSetValue:                           ## @VecSetValue
	.cfi_startproc
## BB#0:
	pushq	%rbp
Ltmp10:
	.cfi_def_cfa_offset 16
Ltmp11:
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
Ltmp12:
	.cfi_def_cfa_register %rbp
	subq	$32, %rsp
	movl	$1, %eax
	leaq	-12(%rbp), %rcx
	leaq	-24(%rbp), %r8
	movq	%rdi, -8(%rbp)
	movl	%esi, -12(%rbp)
	movsd	%xmm0, -24(%rbp)
	movl	%edx, -28(%rbp)
	movq	-8(%rbp), %rdi
	movl	-28(%rbp), %edx
	movl	%eax, %esi
	movl	%edx, -32(%rbp)         ## 4-byte Spill
	movq	%rcx, %rdx
	movq	%r8, %rcx
	movl	-32(%rbp), %r8d         ## 4-byte Reload
	callq	_VecSetValues
	addq	$32, %rsp
	popq	%rbp
	retq
	.cfi_endproc

	.align	4, 0x90
_MatSetValue:                           ## @MatSetValue
	.cfi_startproc
## BB#0:
	pushq	%rbp
Ltmp13:
	.cfi_def_cfa_offset 16
Ltmp14:
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
Ltmp15:
	.cfi_def_cfa_register %rbp
	subq	$48, %rsp
	movl	$1, %eax
	leaq	-12(%rbp), %r8
	leaq	-16(%rbp), %r9
	leaq	-24(%rbp), %r10
	movq	%rdi, -8(%rbp)
	movl	%esi, -12(%rbp)
	movl	%edx, -16(%rbp)
	movsd	%xmm0, -24(%rbp)
	movl	%ecx, -28(%rbp)
	movq	-8(%rbp), %rdi
	movl	-28(%rbp), %ecx
	movl	%eax, %esi
	movq	%r8, %rdx
	movl	%ecx, -32(%rbp)         ## 4-byte Spill
	movl	%eax, %ecx
	movq	%r9, %r8
	movq	%r10, %r9
	movl	-32(%rbp), %eax         ## 4-byte Reload
	movl	%eax, (%rsp)
	callq	_MatSetValues
	addq	$48, %rsp
	popq	%rbp
	retq
	.cfi_endproc

	.section	__TEXT,__cstring,cstring_literals
L_.str:                                 ## @.str
	.asciz	"ConstVars:\n"

L_.str.1:                               ## @.str.1
	.asciz	"%f,%f,%f\n"

L_.str.2:                               ## @.str.2
	.asciz	"%f,%f\n"

L_.str.3:                               ## @.str.3
	.asciz	"%d,%f,%f\n"

L_.str.4:                               ## @.str.4
	.asciz	"Dcs: Ion %d, Comp %d "

L_.str.5:                               ## @.str.5
	.asciz	"Dcs x: %f, Dcs y: %f\n"

L_.str.6:                               ## @.str.6
	.asciz	"Dcb: Ion %d, Comp %d "

L_.str.7:                               ## @.str.7
	.asciz	"Dcb x: %f, Dcb y: %f\n"

L_.str.8:                               ## @.str.8
	.asciz	"Ion: %d, Comp %d, C: %f\n"

L_.str.9:                               ## @.str.9
	.asciz	"Comp %d, Phi: %f\n"

L_.str.10:                              ## @.str.10
	.asciz	"Gvars:\n"

L_.str.11:                              ## @.str.11
	.asciz	"NaT :%f,%f,%f*1e-6\n"

L_.str.12:                              ## @.str.12
	.asciz	"NaP :%f,%f,%f\n"

L_.str.13:                              ## @.str.13
	.asciz	"KDR :%f,%f\n"

L_.str.14:                              ## @.str.14
	.asciz	"KA :%f,%f,%f\n"

L_.str.15:                              ## @.str.15
	.asciz	"Ion: %d, Comp %d\n"

L_.str.16:                              ## @.str.16
	.asciz	"Flux*1e6: %f, dfdci: %f, dfdce: %f, dfdphim: %f\n"

L_.str.17:                              ## @.str.17
	.asciz	"Comp: %d\n"

L_.str.18:                              ## @.str.18
	.asciz	"wFlux: %f,%f,%f\n"

L_.str.19:                              ## @.str.19
	.asciz	"Netwon Iteration did not converge! Stopping...\n"

	.section	__TEXT,__const
	.align	2                       ## @z
_z:
	.long	1                       ## 0x1
	.long	1                       ## 0x1
	.long	4294967295              ## 0xffffffff

	.section	__TEXT,__cstring,cstring_literals
L___func__.calc_residual:               ## @__func__.calc_residual
	.asciz	"calc_residual"

L_.str.20:                              ## @.str.20
	.asciz	"update_solution.c"

L_.str.21:                              ## @.str.21
	.asciz	" "

	.section	__TEXT,__const
	.align	4                       ## @cbath
_cbath:
	.quad	4594212051873190380     ## double 0.14000000000000001
	.quad	4569986288757639007     ## double 0.0033999999999999998
	.quad	4593311331947716280     ## double 0.12

	.section	__TEXT,__cstring,cstring_literals
L_.str.22:                              ## @.str.22
	.asciz	"%f\n"

L_.str.23:                              ## @.str.23
	.asciz	"x:%d,y:%d,ion:%d, %f\n"

	.section	__TEXT,__const
	.align	4                       ## @cm
_cm:
	.quad	4518870892857651536     ## double 1.3264676206595581E-6
	.quad	4518870892857651536     ## double 1.3264676206595581E-6

	.section	__TEXT,__cstring,cstring_literals
L___func__.calc_jacobian:               ## @__func__.calc_jacobian
	.asciz	"calc_jacobian"

L_.str.24:                              ## @.str.24
	.asciz	"Nz: %d, Ind: %d\n"

L_.str.25:                              ## @.str.25
	.asciz	"mat.output"


.subsections_via_symbols
