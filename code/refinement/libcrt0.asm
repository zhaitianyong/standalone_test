	.text	64
	.globl	__crtLibMain
	.extern	_cinit
__crtLibMain:
	cmpl	$1,%edx
	je	_$L1
	movl	$1,%eax
	ret
_$L1:
	callq	_cinit
	movl	$1,%eax
	retq
	.section .rdata
	.long	0,0,0,0,0,0,0,0,0,0
	.long	0,0,0,0,0,0,0,0,0,0
	.long	0,0,0,0,0,0,0,0,0,0
	.long	0,0,0
