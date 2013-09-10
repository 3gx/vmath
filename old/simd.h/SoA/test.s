	.const
__ZStL19piecewise_construct:
	.space 1
	.static_data
__ZStL8__ioinit:
	.space	1
	.const
__ZStL13allocator_arg:
	.space 1
__ZStL6ignore:
	.space 1
	.cstring
LC0:
	.ascii "0\0"
LC1:
	.ascii "test.cpp\0"
	.text
	.align 1,0x90
__ZZ4mainENKUliE_clEi:
LFB3064:
	pushq	%rbp
LCFI0:
	movq	%rsp, %rbp
LCFI1:
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movl	%esi, -12(%rbp)
	movl	-12(%rbp), %eax
	cmpl	$1, %eax
	je	L4
	cmpl	$1, %eax
	jg	L7
	testl	%eax, %eax
	je	L3
	jmp	L2
L7:
	cmpl	$2, %eax
	je	L5
	cmpl	$3, %eax
	je	L6
	jmp	L2
L3:
# 60 "test.cpp" 1
	#case0
# 0 "" 2
	jmp	L1
L4:
# 61 "test.cpp" 1
	#case1
# 0 "" 2
	jmp	L1
L5:
# 62 "test.cpp" 1
	#case2
# 0 "" 2
	jmp	L1
L6:
# 63 "test.cpp" 1
	#case3
# 0 "" 2
	nop
	jmp	L1
L2:
	leaq	LC0(%rip), %rcx
	movl	$64, %edx
	leaq	LC1(%rip), %rsi
	leaq	__ZZZ4mainENKUliE_clEiE8__func__(%rip), %rdi
	call	___assert_rtn
L1:
	leave
LCFI2:
	ret
LFE3064:
	.globl _main
_main:
LFB3063:
	pushq	%rbp
LCFI3:
	movq	%rsp, %rbp
LCFI4:
	subq	$48, %rsp
	movl	%edi, -20(%rbp)
	movq	%rsi, -32(%rbp)
	movl	-20(%rbp), %eax
	movl	%eax, -4(%rbp)
	movb	%dl, (%rsp)
LEHB0:
	call	__ZN10static_forILi0ELi4EE4evalIZ4mainEUliE_EEvT_
LEHE0:
LEHB1:
	call	__Z5printILi10EEvv
	movl	$0, %eax
	jmp	L13
L12:
	movq	%rax, %rdi
	call	__Unwind_Resume
LEHE1:
L13:
	leave
LCFI5:
	ret
LFE3063:
	.section __DATA,__gcc_except_tab
GCC_except_table0:
LLSDA3063:
	.byte	0xff
	.byte	0xff
	.byte	0x3
	.byte	0x1a
	.set L$set$0,LEHB0-LFB3063
	.long L$set$0
	.set L$set$1,LEHE0-LEHB0
	.long L$set$1
	.set L$set$2,L12-LFB3063
	.long L$set$2
	.byte	0
	.set L$set$3,LEHB1-LFB3063
	.long L$set$3
	.set L$set$4,LEHE1-LEHB1
	.long L$set$4
	.long	0
	.byte	0
	.text
__ZN10static_forILi0ELi4EE4evalIZ4mainEUliE_EEvT_:
LFB3130:
	pushq	%rbp
LCFI6:
	movq	%rsp, %rbp
LCFI7:
	pushq	%rbx
	subq	$24, %rsp
LCFI8:
	movl	$0, %esi
	leaq	16(%rbp), %rdi
	call	__ZZ4mainENKUliE_clEi
	movb	%bl, (%rsp)
	call	__ZN10static_forILi1ELi4EE4evalIZ4mainEUliE_EEvT_
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
LCFI9:
	ret
LFE3130:
	.cstring
LC2:
	.ascii " number= %d \12\0"
	.section __TEXT,__textcoal_nt,coalesced,pure_instructions
	.globl __Z5printILi10EEvv
	.weak_definition __Z5printILi10EEvv
__Z5printILi10EEvv:
LFB3131:
	pushq	%rbp
LCFI10:
	movq	%rsp, %rbp
LCFI11:
	movq	___stderrp@GOTPCREL(%rip), %rax
	movq	(%rax), %rax
	movl	$10, %edx
	leaq	LC2(%rip), %rsi
	movq	%rax, %rdi
	movl	$0, %eax
	call	_fprintf
	popq	%rbp
LCFI12:
	ret
LFE3131:
	.text
__ZN10static_forILi1ELi4EE4evalIZ4mainEUliE_EEvT_:
LFB3188:
	pushq	%rbp
LCFI13:
	movq	%rsp, %rbp
LCFI14:
	pushq	%rbx
	subq	$24, %rsp
LCFI15:
	movl	$1, %esi
	leaq	16(%rbp), %rdi
	call	__ZZ4mainENKUliE_clEi
	movb	%bl, (%rsp)
	call	__ZN10static_forILi2ELi4EE4evalIZ4mainEUliE_EEvT_
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
LCFI16:
	ret
LFE3188:
__ZN10static_forILi2ELi4EE4evalIZ4mainEUliE_EEvT_:
LFB3226:
	pushq	%rbp
LCFI17:
	movq	%rsp, %rbp
LCFI18:
	pushq	%rbx
	subq	$24, %rsp
LCFI19:
	movl	$2, %esi
	leaq	16(%rbp), %rdi
	call	__ZZ4mainENKUliE_clEi
	movb	%bl, (%rsp)
	call	__ZN10static_forILi3ELi4EE4evalIZ4mainEUliE_EEvT_
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
LCFI20:
	ret
LFE3226:
__ZN10static_forILi3ELi4EE4evalIZ4mainEUliE_EEvT_:
LFB3237:
	pushq	%rbp
LCFI21:
	movq	%rsp, %rbp
LCFI22:
	pushq	%rbx
	subq	$24, %rsp
LCFI23:
	movl	$3, %esi
	leaq	16(%rbp), %rdi
	call	__ZZ4mainENKUliE_clEi
	movb	%bl, (%rsp)
	call	__ZN10static_forILi4ELi4EE4evalIZ4mainEUliE_EEvT_
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
LCFI24:
	ret
LFE3237:
__ZN10static_forILi4ELi4EE4evalIZ4mainEUliE_EEvT_:
LFB3243:
	pushq	%rbp
LCFI25:
	movq	%rsp, %rbp
LCFI26:
	popq	%rbp
LCFI27:
	ret
LFE3243:
__Z41__static_initialization_and_destruction_0ii:
LFB3262:
	pushq	%rbp
LCFI28:
	movq	%rsp, %rbp
LCFI29:
	subq	$16, %rsp
	movl	%edi, -4(%rbp)
	movl	%esi, -8(%rbp)
	cmpl	$1, -4(%rbp)
	jne	L20
	cmpl	$65535, -8(%rbp)
	jne	L20
	leaq	__ZStL8__ioinit(%rip), %rdi
	call	__ZNSt8ios_base4InitC1Ev
	leaq	___dso_handle(%rip), %rdx
	leaq	__ZStL8__ioinit(%rip), %rsi
	movq	__ZNSt8ios_base4InitD1Ev@GOTPCREL(%rip), %rax
	movq	%rax, %rdi
	call	___cxa_atexit
L20:
	leave
LCFI30:
	ret
LFE3262:
__GLOBAL__sub_I_test.cpp:
LFB3263:
	pushq	%rbp
LCFI31:
	movq	%rsp, %rbp
LCFI32:
	movl	$65535, %esi
	movl	$1, %edi
	call	__Z41__static_initialization_and_destruction_0ii
	popq	%rbp
LCFI33:
	ret
LFE3263:
	.mod_init_func
	.align 3
	.quad	__GLOBAL__sub_I_test.cpp
	.const
__ZZZ4mainENKUliE_clEiE8__func__:
	.ascii "operator()\0"
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
	.quad	LFB3064-.
	.set L$set$7,LFE3064-LFB3064
	.quad L$set$7
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$8,LCFI0-LFB3064
	.long L$set$8
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$9,LCFI1-LCFI0
	.long L$set$9
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$10,LCFI2-LCFI1
	.long L$set$10
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE1:
LSFDE3:
	.set L$set$11,LEFDE3-LASFDE3
	.long L$set$11
LASFDE3:
	.long	LASFDE3-EH_frame1
	.quad	LFB3063-.
	.set L$set$12,LFE3063-LFB3063
	.quad L$set$12
	.byte	0x8
	.quad	LLSDA3063-.
	.byte	0x4
	.set L$set$13,LCFI3-LFB3063
	.long L$set$13
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$14,LCFI4-LCFI3
	.long L$set$14
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$15,LCFI5-LCFI4
	.long L$set$15
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE3:
LSFDE5:
	.set L$set$16,LEFDE5-LASFDE5
	.long L$set$16
LASFDE5:
	.long	LASFDE5-EH_frame1
	.quad	LFB3130-.
	.set L$set$17,LFE3130-LFB3130
	.quad L$set$17
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$18,LCFI6-LFB3130
	.long L$set$18
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$19,LCFI7-LCFI6
	.long L$set$19
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$20,LCFI8-LCFI7
	.long L$set$20
	.byte	0x83
	.byte	0x3
	.byte	0x4
	.set L$set$21,LCFI9-LCFI8
	.long L$set$21
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE5:
LSFDE7:
	.set L$set$22,LEFDE7-LASFDE7
	.long L$set$22
LASFDE7:
	.long	LASFDE7-EH_frame1
	.quad	LFB3131-.
	.set L$set$23,LFE3131-LFB3131
	.quad L$set$23
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$24,LCFI10-LFB3131
	.long L$set$24
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$25,LCFI11-LCFI10
	.long L$set$25
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$26,LCFI12-LCFI11
	.long L$set$26
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE7:
LSFDE9:
	.set L$set$27,LEFDE9-LASFDE9
	.long L$set$27
LASFDE9:
	.long	LASFDE9-EH_frame1
	.quad	LFB3188-.
	.set L$set$28,LFE3188-LFB3188
	.quad L$set$28
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$29,LCFI13-LFB3188
	.long L$set$29
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$30,LCFI14-LCFI13
	.long L$set$30
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$31,LCFI15-LCFI14
	.long L$set$31
	.byte	0x83
	.byte	0x3
	.byte	0x4
	.set L$set$32,LCFI16-LCFI15
	.long L$set$32
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE9:
LSFDE11:
	.set L$set$33,LEFDE11-LASFDE11
	.long L$set$33
LASFDE11:
	.long	LASFDE11-EH_frame1
	.quad	LFB3226-.
	.set L$set$34,LFE3226-LFB3226
	.quad L$set$34
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$35,LCFI17-LFB3226
	.long L$set$35
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$36,LCFI18-LCFI17
	.long L$set$36
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$37,LCFI19-LCFI18
	.long L$set$37
	.byte	0x83
	.byte	0x3
	.byte	0x4
	.set L$set$38,LCFI20-LCFI19
	.long L$set$38
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE11:
LSFDE13:
	.set L$set$39,LEFDE13-LASFDE13
	.long L$set$39
LASFDE13:
	.long	LASFDE13-EH_frame1
	.quad	LFB3237-.
	.set L$set$40,LFE3237-LFB3237
	.quad L$set$40
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$41,LCFI21-LFB3237
	.long L$set$41
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$42,LCFI22-LCFI21
	.long L$set$42
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$43,LCFI23-LCFI22
	.long L$set$43
	.byte	0x83
	.byte	0x3
	.byte	0x4
	.set L$set$44,LCFI24-LCFI23
	.long L$set$44
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE13:
LSFDE15:
	.set L$set$45,LEFDE15-LASFDE15
	.long L$set$45
LASFDE15:
	.long	LASFDE15-EH_frame1
	.quad	LFB3243-.
	.set L$set$46,LFE3243-LFB3243
	.quad L$set$46
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$47,LCFI25-LFB3243
	.long L$set$47
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$48,LCFI26-LCFI25
	.long L$set$48
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$49,LCFI27-LCFI26
	.long L$set$49
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE15:
LSFDE17:
	.set L$set$50,LEFDE17-LASFDE17
	.long L$set$50
LASFDE17:
	.long	LASFDE17-EH_frame1
	.quad	LFB3262-.
	.set L$set$51,LFE3262-LFB3262
	.quad L$set$51
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$52,LCFI28-LFB3262
	.long L$set$52
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$53,LCFI29-LCFI28
	.long L$set$53
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$54,LCFI30-LCFI29
	.long L$set$54
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE17:
LSFDE19:
	.set L$set$55,LEFDE19-LASFDE19
	.long L$set$55
LASFDE19:
	.long	LASFDE19-EH_frame1
	.quad	LFB3263-.
	.set L$set$56,LFE3263-LFB3263
	.quad L$set$56
	.byte	0x8
	.quad	0
	.byte	0x4
	.set L$set$57,LCFI31-LFB3263
	.long L$set$57
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$58,LCFI32-LCFI31
	.long L$set$58
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$59,LCFI33-LCFI32
	.long L$set$59
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE19:
	.constructor
	.destructor
	.align 1
	.subsections_via_symbols
