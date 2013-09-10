	.section __TEXT,__textcoal_nt,coalesced,pure_instructions
	.align 1
	.align 4
	.globl __ZN12ParticleDataIdED1Ev
	.weak_definition __ZN12ParticleDataIdED1Ev
__ZN12ParticleDataIdED1Ev:
LFB1754:
	pushq	%rbx
LCFI0:
	movq	%rdi, %rbx
	movq	144(%rdi), %rdi
	testq	%rdi, %rdi
	je	L2
	call	__ZdlPv
L2:
	movq	120(%rbx), %rdi
	testq	%rdi, %rdi
	je	L3
	call	__ZdlPv
L3:
	movq	96(%rbx), %rdi
	testq	%rdi, %rdi
	je	L4
	call	__ZdlPv
L4:
	movq	72(%rbx), %rdi
	testq	%rdi, %rdi
	je	L5
	call	__ZdlPv
L5:
	movq	48(%rbx), %rdi
	testq	%rdi, %rdi
	je	L6
	call	__ZdlPv
L6:
	movq	24(%rbx), %rdi
	testq	%rdi, %rdi
	je	L7
	call	__ZdlPv
L7:
	movq	(%rbx), %rdi
	testq	%rdi, %rdi
	je	L1
	popq	%rbx
LCFI1:
	jmp	__ZdlPv
	.align 4
L1:
LCFI2:
	popq	%rbx
LCFI3:
	ret
LFE1754:
	.cstring
LC1:
	.ascii "vector::_M_default_append\0"
	.section __TEXT,__textcoal_nt,coalesced,pure_instructions
	.align 1
	.align 4
	.globl __ZNSt6vectorIdSaIdEE17_M_default_appendEm
	.weak_definition __ZNSt6vectorIdSaIdEE17_M_default_appendEm
__ZNSt6vectorIdSaIdEE17_M_default_appendEm:
LFB1931:
	pushq	%r14
LCFI4:
	testq	%rsi, %rsi
	pushq	%r13
LCFI5:
	pushq	%r12
LCFI6:
	pushq	%rbp
LCFI7:
	pushq	%rbx
LCFI8:
	movq	%rsi, %rbx
	je	L29
	movq	8(%rdi), %r10
	movq	%rdi, %rbp
	movq	16(%rdi), %rax
	subq	%r10, %rax
	sarq	$3, %rax
	cmpq	%rax, %rsi
	ja	L31
	movq	%r10, %rsi
	salq	$60, %rsi
	shrq	$63, %rsi
	cmpq	%rbx, %rsi
	cmova	%rbx, %rsi
	cmpq	$6, %rbx
	ja	L82
	movq	%rbx, %rsi
L32:
	movq	%rbx, %rdx
	movq	%r10, %rax
	xorl	%ecx, %ecx
	.align 4
L34:
	addq	$1, %rcx
	movq	$0, (%rax)
	subq	$1, %rdx
	addq	$8, %rax
	cmpq	%rcx, %rsi
	ja	L34
	cmpq	%rsi, %rbx
	je	L35
L33:
	movq	%rbx, %r11
	subq	%rsi, %r11
	movq	%r11, %r9
	shrq	%r9
	movq	%r9, %rdi
	addq	%rdi, %rdi
	je	L69
	leaq	(%r10,%rsi,8), %r8
	xorl	%ecx, %ecx
	xorpd	%xmm0, %xmm0
	.align 4
L37:
	addq	$1, %rcx
	movapd	%xmm0, (%r8)
	addq	$16, %r8
	cmpq	%r9, %rcx
	jb	L37
	leaq	(%rax,%rdi,8), %rax
	subq	%rdi, %rdx
	cmpq	%rdi, %r11
	je	L35
	.align 4
L69:
	movq	$0, (%rax)
	addq	$8, %rax
	subq	$1, %rdx
	jne	L69
L35:
	leaq	(%r10,%rbx,8), %rax
	movq	%rax, 8(%rbp)
L29:
	popq	%rbx
LCFI9:
	popq	%rbp
LCFI10:
	popq	%r12
LCFI11:
	popq	%r13
LCFI12:
	popq	%r14
LCFI13:
	ret
	.align 4
L31:
LCFI14:
	movabsq	$2305843009213693951, %rsi
	movq	(%rdi), %rdi
	movq	%rsi, %rdx
	subq	%rdi, %r10
	sarq	$3, %r10
	subq	%r10, %rdx
	movq	%r10, %rax
	movq	%r10, %rcx
	cmpq	%rdx, %rbx
	ja	L83
	cmpq	%rbx, %r10
	movq	%rbx, %rdx
	movq	$-8, %r12
	cmovae	%r10, %rdx
	addq	%rdx, %r10
	jae	L84
L41:
	movq	%r12, %rdi
	call	__Znwm
	movq	0(%rbp), %rdi
	movq	%rax, %r13
	movq	8(%rbp), %rax
	subq	%rdi, %rax
	sarq	$3, %rax
	movq	%rax, %rcx
L42:
	leaq	0(,%rcx,8), %r14
	testq	%rax, %rax
	je	L44
	movq	%rdi, %rsi
	movq	%r14, %rdx
	movq	%r13, %rdi
	call	_memmove
	movq	0(%rbp), %rdi
L44:
	leaq	0(%r13,%r14), %rax
	movq	%rax, %r9
	salq	$60, %r9
	shrq	$63, %r9
	cmpq	%rbx, %r9
	cmova	%rbx, %r9
	cmpq	$6, %rbx
	ja	L85
	movq	%rbx, %r9
L45:
	movq	%rbx, %rcx
	movq	%rax, %rdx
	xorl	%r8d, %r8d
	.align 4
L47:
	addq	$1, %r8
	movq	$0, (%rdx)
	subq	$1, %rcx
	addq	$8, %rdx
	cmpq	%r8, %r9
	ja	L47
	cmpq	%r9, %rbx
	je	L48
L46:
	movq	%rbx, %r11
	subq	%r9, %r11
	movq	%r11, %r10
	shrq	%r10
	movq	%r10, %rsi
	addq	%rsi, %rsi
	je	L70
	leaq	(%rax,%r9,8), %r9
	xorl	%r8d, %r8d
	xorpd	%xmm0, %xmm0
	.align 4
L50:
	addq	$1, %r8
	movapd	%xmm0, (%r9)
	addq	$16, %r9
	cmpq	%r10, %r8
	jb	L50
	leaq	(%rdx,%rsi,8), %rdx
	subq	%rsi, %rcx
	cmpq	%rsi, %r11
	je	L48
	.align 4
L70:
	movq	$0, (%rdx)
	addq	$8, %rdx
	subq	$1, %rcx
	jne	L70
L48:
	leaq	(%rax,%rbx,8), %rbx
	testq	%rdi, %rdi
	je	L53
	call	__ZdlPv
L53:
	addq	%r13, %r12
	movq	%r13, 0(%rbp)
	movq	%rbx, 8(%rbp)
	movq	%r12, 16(%rbp)
	popq	%rbx
LCFI15:
	popq	%rbp
LCFI16:
	popq	%r12
LCFI17:
	popq	%r13
LCFI18:
	popq	%r14
LCFI19:
	ret
	.align 4
L84:
LCFI20:
	cmpq	%rsi, %r10
	ja	L41
	testq	%r10, %r10
	jne	L86
	xorl	%r12d, %r12d
	xorl	%r13d, %r13d
	jmp	L42
L85:
	testq	%r9, %r9
	jne	L45
	movq	%rbx, %rcx
	movq	%rax, %rdx
	jmp	L46
L86:
	leaq	0(,%r10,8), %r12
	jmp	L41
L83:
	leaq	LC1(%rip), %rdi
	call	__ZSt20__throw_length_errorPKc
L82:
	testq	%rsi, %rsi
	jne	L32
	movq	%rbx, %rdx
	movq	%r10, %rax
	jmp	L33
LFE1931:
	.align 1
	.align 4
	.globl __ZN12ParticleDataIdE6resizeEm
	.weak_definition __ZN12ParticleDataIdE6resizeEm
__ZN12ParticleDataIdE6resizeEm:
LFB1822:
	pushq	%rbp
LCFI21:
	movq	%rsi, %rbp
	pushq	%rbx
LCFI22:
	movq	%rdi, %rbx
	subq	$8, %rsp
LCFI23:
	movq	(%rdi), %rdx
	movq	8(%rdi), %rax
	subq	%rdx, %rax
	sarq	$3, %rax
	cmpq	%rax, %rsi
	ja	L103
	jae	L89
	leaq	(%rdx,%rsi,8), %rax
	movq	%rax, 8(%rdi)
L89:
	movq	24(%rbx), %rdx
	movq	32(%rbx), %rax
	subq	%rdx, %rax
	sarq	$3, %rax
	cmpq	%rax, %rbp
	ja	L104
	jae	L91
	leaq	(%rdx,%rbp,8), %rax
	movq	%rax, 32(%rbx)
L91:
	movq	48(%rbx), %rdx
	movq	56(%rbx), %rax
	subq	%rdx, %rax
	sarq	$3, %rax
	cmpq	%rax, %rbp
	ja	L105
	jae	L93
	leaq	(%rdx,%rbp,8), %rax
	movq	%rax, 56(%rbx)
L93:
	movq	72(%rbx), %rdx
	movq	80(%rbx), %rax
	subq	%rdx, %rax
	sarq	$3, %rax
	cmpq	%rax, %rbp
	ja	L106
	jae	L95
	leaq	(%rdx,%rbp,8), %rax
	movq	%rax, 80(%rbx)
L95:
	movq	96(%rbx), %rdx
	movq	104(%rbx), %rax
	subq	%rdx, %rax
	sarq	$3, %rax
	cmpq	%rax, %rbp
	ja	L107
	jae	L97
	leaq	(%rdx,%rbp,8), %rax
	movq	%rax, 104(%rbx)
L97:
	movq	120(%rbx), %rdx
	movq	128(%rbx), %rax
	subq	%rdx, %rax
	sarq	$3, %rax
	cmpq	%rax, %rbp
	ja	L108
	jae	L99
	leaq	(%rdx,%rbp,8), %rax
	movq	%rax, 128(%rbx)
L99:
	movq	144(%rbx), %rdx
	movq	152(%rbx), %rax
	subq	%rdx, %rax
	sarq	$3, %rax
	cmpq	%rax, %rbp
	ja	L109
	jae	L87
	leaq	(%rdx,%rbp,8), %rax
	movq	%rax, 152(%rbx)
L87:
	addq	$8, %rsp
LCFI24:
	popq	%rbx
LCFI25:
	popq	%rbp
LCFI26:
	ret
	.align 4
L103:
LCFI27:
	subq	%rax, %rsi
	call	__ZNSt6vectorIdSaIdEE17_M_default_appendEm
	jmp	L89
	.align 4
L105:
	leaq	48(%rbx), %rdi
	movq	%rbp, %rsi
	subq	%rax, %rsi
	call	__ZNSt6vectorIdSaIdEE17_M_default_appendEm
	jmp	L93
	.align 4
L104:
	leaq	24(%rbx), %rdi
	movq	%rbp, %rsi
	subq	%rax, %rsi
	call	__ZNSt6vectorIdSaIdEE17_M_default_appendEm
	jmp	L91
	.align 4
L109:
	leaq	144(%rbx), %rdi
	addq	$8, %rsp
LCFI28:
	movq	%rbp, %rsi
	popq	%rbx
LCFI29:
	subq	%rax, %rsi
	popq	%rbp
LCFI30:
	jmp	__ZNSt6vectorIdSaIdEE17_M_default_appendEm
	.align 4
L108:
LCFI31:
	leaq	120(%rbx), %rdi
	movq	%rbp, %rsi
	subq	%rax, %rsi
	call	__ZNSt6vectorIdSaIdEE17_M_default_appendEm
	jmp	L99
	.align 4
L107:
	leaq	96(%rbx), %rdi
	movq	%rbp, %rsi
	subq	%rax, %rsi
	call	__ZNSt6vectorIdSaIdEE17_M_default_appendEm
	jmp	L97
	.align 4
L106:
	leaq	72(%rbx), %rdi
	movq	%rbp, %rsi
	subq	%rax, %rsi
	call	__ZNSt6vectorIdSaIdEE17_M_default_appendEm
	jmp	L95
LFE1822:
	.cstring
LC3:
	.ascii "mTotal= %g\12\0"
LC6:
	.ascii "myvec= %g %g %g \12\0"
LC8:
	.ascii ">>mTotal= %g\12\0"
	.section __TEXT,__text_startup,regular,pure_instructions
	.align 4
	.globl _main
_main:
LFB1747:
	pushq	%r12
LCFI32:
	movl	$1024, %esi
	pushq	%rbp
LCFI33:
	movl	%edi, %ebp
	pushq	%rbx
LCFI34:
	subq	$224, %rsp
LCFI35:
	movq	$0, 48(%rsp)
	leaq	48(%rsp), %rdi
	movq	$0, 56(%rsp)
	movq	$0, 64(%rsp)
	movq	$0, 72(%rsp)
	movq	$0, 80(%rsp)
	movq	$0, 88(%rsp)
	movq	$0, 96(%rsp)
	movq	$0, 104(%rsp)
	movq	$0, 112(%rsp)
	movq	$0, 120(%rsp)
	movq	$0, 128(%rsp)
	movq	$0, 136(%rsp)
	movq	$0, 144(%rsp)
	movq	$0, 152(%rsp)
	movq	$0, 160(%rsp)
	movq	$0, 168(%rsp)
	movq	$0, 176(%rsp)
	movq	$0, 184(%rsp)
	movq	$0, 192(%rsp)
	movq	$0, 200(%rsp)
	movq	$0, 208(%rsp)
LEHB0:
	call	__ZN12ParticleDataIdE6resizeEm
	xorl	%ebx, %ebx
	.align 4
L111:
	call	_drand48
	movsd	%xmm0, 8(%rsp)
	call	_drand48
	movsd	%xmm0, 16(%rsp)
	call	_drand48
	movsd	%xmm0, 24(%rsp)
	call	_drand48
	movsd	%xmm0, 32(%rsp)
	call	_drand48
	movsd	%xmm0, 40(%rsp)
	call	_drand48
	movq	48(%rsp), %rdx
	movq	72(%rsp), %r9
	movsd	24(%rsp), %xmm1
	movq	120(%rsp), %rdi
	movq	96(%rsp), %r8
	movq	192(%rsp), %rax
	addq	%rbx, %r9
	movq	168(%rsp), %rcx
	movq	144(%rsp), %rsi
	movsd	%xmm1, (%rdx,%rbx)
	addq	%rbx, %rdi
	movsd	16(%rsp), %xmm1
	addq	%rbx, %r8
	addq	%rbx, %rax
	movsd	%xmm1, (%r9)
	movsd	8(%rsp), %xmm1
	addq	%rbx, %rcx
	addq	%rbx, %rsi
	addq	$8, %rbx
	movsd	%xmm1, (%r8)
	cmpq	$8192, %rbx
	movsd	%xmm0, (%rdi)
	movsd	40(%rsp), %xmm0
	movsd	%xmm0, (%rsi)
	movsd	32(%rsp), %xmm0
	movsd	%xmm0, (%rcx)
	movabsq	$4562146422526312448, %rcx
	movq	%rcx, (%rax)
	jne	L111
	xorl	%eax, %eax
	xorl	%ebx, %ebx
	.align 4
L112:
	movd	%rbx, %xmm1
	addsd	(%rdx,%rax), %xmm1
	addq	$8, %rax
	cmpq	$8192, %rax
	movd	%xmm1, %rbx
	jne	L112
	movq	___stderrp@GOTPCREL(%rip), %r12
	movapd	%xmm1, %xmm0
	movl	$1, %eax
	leaq	LC3(%rip), %rsi
	movq	(%r12), %rdi
	call	_fprintf
	cvtsi2sd	%ebp, %xmm0
	movq	(%r12), %rdi
	movl	$3, %eax
	movsd	LC4(%rip), %xmm2
	leaq	LC6(%rip), %rsi
	movsd	LC5(%rip), %xmm1
	movsd	%xmm0, 8(%rsp)
	call	_fprintf
# 58 "nbody_lf.cpp" 1
	#begin-A
# 0 "" 2
	movsd	8(%rsp), %xmm0
	mulsd	%xmm0, %xmm0
	addsd	8(%rsp), %xmm0
# 75 "nbody_lf.cpp" 1
	#end-A
# 0 "" 2
# 76 "nbody_lf.cpp" 1
	#begin-B
# 0 "" 2
	movsd	LC7(%rip), %xmm1
	mulsd	%xmm0, %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, 8(%rsp)
# 81 "nbody_lf.cpp" 1
	#end-B
# 0 "" 2
	movq	(%r12), %rdi
	movd	%rbx, %xmm0
	movl	$1, %eax
	leaq	LC8(%rip), %rsi
	call	_fprintf
	movq	(%r12), %rdi
	leaq	LC6(%rip), %rsi
	movl	$3, %eax
	movsd	LC9(%rip), %xmm2
	movsd	LC10(%rip), %xmm1
	movsd	8(%rsp), %xmm0
	call	_fprintf
LEHE0:
	leaq	48(%rsp), %rdi
	call	__ZN12ParticleDataIdED1Ev
	addq	$224, %rsp
LCFI36:
	xorl	%eax, %eax
	popq	%rbx
LCFI37:
	popq	%rbp
LCFI38:
	popq	%r12
LCFI39:
	ret
L114:
LCFI40:
	leaq	48(%rsp), %rdi
	movq	%rax, %rbx
	call	__ZN12ParticleDataIdED1Ev
	movq	%rbx, %rdi
LEHB1:
	call	__Unwind_Resume
LEHE1:
LFE1747:
	.section __DATA,__gcc_except_tab
GCC_except_table0:
LLSDA1747:
	.byte	0xff
	.byte	0xff
	.byte	0x3
	.byte	0x1a
	.set L$set$0,LEHB0-LFB1747
	.long L$set$0
	.set L$set$1,LEHE0-LEHB0
	.long L$set$1
	.set L$set$2,L114-LFB1747
	.long L$set$2
	.byte	0
	.set L$set$3,LEHB1-LFB1747
	.long L$set$3
	.set L$set$4,LEHE1-LEHB1
	.long L$set$4
	.long	0
	.byte	0
	.section __TEXT,__text_startup,regular,pure_instructions
	.align 4
__GLOBAL__sub_I_nbody_lf.cpp:
LFB2002:
	leaq	__ZStL8__ioinit(%rip), %rdi
	subq	$8, %rsp
LCFI41:
	call	__ZNSt8ios_base4InitC1Ev
	movq	__ZNSt8ios_base4InitD1Ev@GOTPCREL(%rip), %rdi
	addq	$8, %rsp
LCFI42:
	leaq	___dso_handle(%rip), %rdx
	leaq	__ZStL8__ioinit(%rip), %rsi
	jmp	___cxa_atexit
LFE2002:
	.mod_init_func
	.align 3
	.quad	__GLOBAL__sub_I_nbody_lf.cpp
	.static_data
__ZStL8__ioinit:
	.space	1
	.literal8
	.align 3
LC4:
	.long	0
	.long	1074266112
	.align 3
LC5:
	.long	0
	.long	1073741824
	.align 3
LC7:
	.long	0
	.long	1075314688
	.align 3
LC9:
	.long	0
	.long	1076363264
	.align 3
LC10:
	.long	0
	.long	1077411840
	.section __TEXT,__eh_frame,coalesced,no_toc+strip_static_syms+live_support
EH_frame1:
	.set L$set$5,LECIE1-LSCIE1
	.long L$set$5
LSCIE1:
	.long	0
	.byte	0x1
	.ascii "zPLR\0"
	.byte	0x1
	.byte	0x78
	.byte	0x10
	.byte	0x7
	.byte	0x9b
	.long	___gxx_personality_v0+4@GOTPCREL
	.byte	0x10
	.byte	0x10
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.byte	0x90
	.byte	0x1
	.align 3
LECIE1:
LSFDE1:
	.set L$set$6,LEFDE1-LASFDE1
	.long L$set$6
LASFDE1:
	.long	LASFDE1-EH_frame1
	.quad	LFB1754-.
	.set L$set$7,LFE1754-LFB1754
	.quad L$set$7
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$8,LCFI0-LFB1754
	.long L$set$8
	.byte	0xe
	.byte	0x10
	.byte	0x83
	.byte	0x2
	.byte	0x4
	.set L$set$9,LCFI1-LCFI0
	.long L$set$9
	.byte	0xa
	.byte	0xe
	.byte	0x8
	.byte	0x4
	.set L$set$10,LCFI2-LCFI1
	.long L$set$10
	.byte	0xb
	.byte	0x4
	.set L$set$11,LCFI3-LCFI2
	.long L$set$11
	.byte	0xe
	.byte	0x8
	.align 3
LEFDE1:
LSFDE3:
	.set L$set$12,LEFDE3-LASFDE3
	.long L$set$12
LASFDE3:
	.long	LASFDE3-EH_frame1
	.quad	LFB1931-.
	.set L$set$13,LFE1931-LFB1931
	.quad L$set$13
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$14,LCFI4-LFB1931
	.long L$set$14
	.byte	0xe
	.byte	0x10
	.byte	0x8e
	.byte	0x2
	.byte	0x4
	.set L$set$15,LCFI5-LCFI4
	.long L$set$15
	.byte	0xe
	.byte	0x18
	.byte	0x8d
	.byte	0x3
	.byte	0x4
	.set L$set$16,LCFI6-LCFI5
	.long L$set$16
	.byte	0xe
	.byte	0x20
	.byte	0x8c
	.byte	0x4
	.byte	0x4
	.set L$set$17,LCFI7-LCFI6
	.long L$set$17
	.byte	0xe
	.byte	0x28
	.byte	0x86
	.byte	0x5
	.byte	0x4
	.set L$set$18,LCFI8-LCFI7
	.long L$set$18
	.byte	0xe
	.byte	0x30
	.byte	0x83
	.byte	0x6
	.byte	0x4
	.set L$set$19,LCFI9-LCFI8
	.long L$set$19
	.byte	0xa
	.byte	0xe
	.byte	0x28
	.byte	0x4
	.set L$set$20,LCFI10-LCFI9
	.long L$set$20
	.byte	0xe
	.byte	0x20
	.byte	0x4
	.set L$set$21,LCFI11-LCFI10
	.long L$set$21
	.byte	0xe
	.byte	0x18
	.byte	0x4
	.set L$set$22,LCFI12-LCFI11
	.long L$set$22
	.byte	0xe
	.byte	0x10
	.byte	0x4
	.set L$set$23,LCFI13-LCFI12
	.long L$set$23
	.byte	0xe
	.byte	0x8
	.byte	0x4
	.set L$set$24,LCFI14-LCFI13
	.long L$set$24
	.byte	0xb
	.byte	0x4
	.set L$set$25,LCFI15-LCFI14
	.long L$set$25
	.byte	0xa
	.byte	0xe
	.byte	0x28
	.byte	0x4
	.set L$set$26,LCFI16-LCFI15
	.long L$set$26
	.byte	0xe
	.byte	0x20
	.byte	0x4
	.set L$set$27,LCFI17-LCFI16
	.long L$set$27
	.byte	0xe
	.byte	0x18
	.byte	0x4
	.set L$set$28,LCFI18-LCFI17
	.long L$set$28
	.byte	0xe
	.byte	0x10
	.byte	0x4
	.set L$set$29,LCFI19-LCFI18
	.long L$set$29
	.byte	0xe
	.byte	0x8
	.byte	0x4
	.set L$set$30,LCFI20-LCFI19
	.long L$set$30
	.byte	0xb
	.align 3
LEFDE3:
LSFDE5:
	.set L$set$31,LEFDE5-LASFDE5
	.long L$set$31
LASFDE5:
	.long	LASFDE5-EH_frame1
	.quad	LFB1822-.
	.set L$set$32,LFE1822-LFB1822
	.quad L$set$32
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$33,LCFI21-LFB1822
	.long L$set$33
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$34,LCFI22-LCFI21
	.long L$set$34
	.byte	0xe
	.byte	0x18
	.byte	0x83
	.byte	0x3
	.byte	0x4
	.set L$set$35,LCFI23-LCFI22
	.long L$set$35
	.byte	0xe
	.byte	0x20
	.byte	0x4
	.set L$set$36,LCFI24-LCFI23
	.long L$set$36
	.byte	0xa
	.byte	0xe
	.byte	0x18
	.byte	0x4
	.set L$set$37,LCFI25-LCFI24
	.long L$set$37
	.byte	0xe
	.byte	0x10
	.byte	0x4
	.set L$set$38,LCFI26-LCFI25
	.long L$set$38
	.byte	0xe
	.byte	0x8
	.byte	0x4
	.set L$set$39,LCFI27-LCFI26
	.long L$set$39
	.byte	0xb
	.byte	0x4
	.set L$set$40,LCFI28-LCFI27
	.long L$set$40
	.byte	0xa
	.byte	0xe
	.byte	0x18
	.byte	0x4
	.set L$set$41,LCFI29-LCFI28
	.long L$set$41
	.byte	0xe
	.byte	0x10
	.byte	0x4
	.set L$set$42,LCFI30-LCFI29
	.long L$set$42
	.byte	0xe
	.byte	0x8
	.byte	0x4
	.set L$set$43,LCFI31-LCFI30
	.long L$set$43
	.byte	0xb
	.align 3
LEFDE5:
LSFDE7:
	.set L$set$44,LEFDE7-LASFDE7
	.long L$set$44
LASFDE7:
	.long	LASFDE7-EH_frame1
	.quad	LFB1747-.
	.set L$set$45,LFE1747-LFB1747
	.quad L$set$45
	.byte	0x8
	.quad	LLSDA1747-.
	.byte	0x4
	.set L$set$46,LCFI32-LFB1747
	.long L$set$46
	.byte	0xe
	.byte	0x10
	.byte	0x8c
	.byte	0x2
	.byte	0x4
	.set L$set$47,LCFI33-LCFI32
	.long L$set$47
	.byte	0xe
	.byte	0x18
	.byte	0x86
	.byte	0x3
	.byte	0x4
	.set L$set$48,LCFI34-LCFI33
	.long L$set$48
	.byte	0xe
	.byte	0x20
	.byte	0x83
	.byte	0x4
	.byte	0x4
	.set L$set$49,LCFI35-LCFI34
	.long L$set$49
	.byte	0xe
	.byte	0x80,0x2
	.byte	0x4
	.set L$set$50,LCFI36-LCFI35
	.long L$set$50
	.byte	0xa
	.byte	0xe
	.byte	0x20
	.byte	0x4
	.set L$set$51,LCFI37-LCFI36
	.long L$set$51
	.byte	0xe
	.byte	0x18
	.byte	0x4
	.set L$set$52,LCFI38-LCFI37
	.long L$set$52
	.byte	0xe
	.byte	0x10
	.byte	0x4
	.set L$set$53,LCFI39-LCFI38
	.long L$set$53
	.byte	0xe
	.byte	0x8
	.byte	0x4
	.set L$set$54,LCFI40-LCFI39
	.long L$set$54
	.byte	0xb
	.align 3
LEFDE7:
LSFDE9:
	.set L$set$55,LEFDE9-LASFDE9
	.long L$set$55
LASFDE9:
	.long	LASFDE9-EH_frame1
	.quad	LFB2002-.
	.set L$set$56,LFE2002-LFB2002
	.quad L$set$56
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$57,LCFI41-LFB2002
	.long L$set$57
	.byte	0xe
	.byte	0x10
	.byte	0x4
	.set L$set$58,LCFI42-LCFI41
	.long L$set$58
	.byte	0xe
	.byte	0x8
	.align 3
LEFDE9:
	.constructor
	.destructor
	.align 1
	.subsections_via_symbols
