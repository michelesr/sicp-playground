(define nil '())
(define (square x) (* x x))
(define (cube x) (* x x x))

(define (average a b)
  (/ (+ a b) 2))

(define (smallest-divisor n) (find-divisor n 2))

(define (find-divisor n test-divisor)
  (cond ((> (square test-divisor) n) n)
        ((divides? test-divisor n) test-divisor)
        (else (find-divisor n (+ test-divisor 1)))))
(define (divides? a b) (= (remainder b a) 0))

(define (prime? n)
(= n (smallest-divisor n)))

(define (gcd a b)
  (if (= b 0)
    a
    (gcd b (remainder a b))))

;; (define (sum term a next b)
;;   (if (> a b)
;;     0
;;     (+ (term a)
;;        (sum term (next a) next b))))

(define (product term a next b)
  (if (> a b)
    1
    (* (term a)
       (product term (next a) next b))))

(define (id x) x)
(define (inc x) (+ 1 x))

;; (define (sum term a next b)
;;   (define (iter a result)
;;     (if (> a b)
;;       result
;;       (iter (next a) (+ result (term a)))))
;;   (iter a 0))
;;
;; (define (product term a next b)
;;   (define (iter a result)
;;     (if (> a b)
;;       result
;;       (iter (next a) (* result (term a)))))
;;   (iter a 1))

(define (accumulate combiner null-value term a next b)
  (define (iter a result)
    (if (> a b)
      result
      (iter (next a) (combiner result (term a)))))
  (iter a null-value))

(define (filtered-accumulate combiner filter null-value term a next b)
  (define (iter a result)
    (if (> a b)
      result
      (iter (next a)
            (combiner
              result
              (if (filter a) (term a) null-value)))))
  (iter a null-value))

;; (define (accumulate combiner null-value term a next b)
;;   (if (> a b)
;;     null-value
;;     (combiner (term a)
;;        (accumulate combiner null-value term (next a) next b))))

(define (product term a next b) (accumulate * 1 term a next b))
(define (sum term a next b) (accumulate + 0 term a next b))

(define (sum-of-prime-squares a b)
  (filtered-accumulate + prime? 0 square a inc b))

(define (product-of-relative-primes n)
  (define (is-relative-prime? i) (= (gcd i n) 1))
  (filtered-accumulate * is-relative-prime? 1 id 1 inc n))

(define (pi n)
  (define (term n)
    (/ (* 4 (square n))
       (- (* 4 (square n)) 1)))
  (* 2 (product term 1.0 inc n)))

(define (integral f a b dx)
  (define (add-dx x)
    (+ x dx))
  (* (sum f (+ a (/ dx 2.0)) add-dx b)
     dx))

; WARNING: n must be even!
(define (simpson-integral f a b n)
  (let ((h (/ (- b a) n)))
    (define (next k) (+ k 2))
    (define (term k)
      (+ (* 4 (f (+ a (* k h))))
         (* 2 (f (+ a (* (+ k 1) h))))))
    (/ (* h (+ (f a) (sum term 1 next (- n 1)))) 3)))


(define (frac n d k)
  (define (recur i)
    (/ (n i)
       (+ (d i) (if (>= i k) 0 (recur (+ i 1))))))
  (recur 1))

(define (n i) 1)
(define (d i)
  (let ((r (remainder i 3)))
    (cond ((= r 1)
           (let ((j (/ (+ i 2) 3)))
             (* 2 j)))
          ((= r 2) 1)
          ((= r 0) 1))))

(define (frac-i n d k)
  (define (iter i res)
    (if (< i 0)
      res
      (iter (- i 1) (/ (n i) (+ res (d i))))))
  (iter k 0))

(define (tan x k)
  (frac-i
   (lambda (i) (if (= (+ i 1) 1) x (* x x -1)))
   (lambda (i) (- (* 2.0 (+ i 1)) 1))
   k))

(define (average-damp f)
(lambda (x) (average x (f x))))

(define dx 0.00001)

(define (deriv g)
  (lambda (x) (/ (- (g (+ x dx)) (g x)) dx)))

(define (newton-transform g)
  (lambda (x) (- x (/ (g x) ((deriv g) x)))))
(define (newtons-method g guess)
  (fixed-point (newton-transform g) guess))

(define (fixed-point-of-transform g transform guess)
  (fixed-point (transform g) guess))

(define (cubic a b c)
  (lambda (x) (+ (* x x x) (* a x x) (* b x) c)))

(define (double f) (lambda (x) (f (f x))))
(define (compose f g) (lambda (x) (f (g x))))

(define (repeated f n)
  (if (= n 1) f (compose f (repeated f (- n 1)))))

(define (smooth f)
  (lambda (x)
    (/ (+
         (f (- x dx))
         (f x)
         (f (+ x dx)))
       3)))

(define (n-smooth f n)
  ((repeated smooth n) f))

;; (define (sqrt x)
;;   (sqrt-iter 1.0 x))
;;
;; (define (sqrt-iter guess x)
;;   (if (good-enough? guess x)
;;     guess
;;     (sqrt-iter (improve guess x) x)))
;; (define (good-enough? guess x)
;;   (< (abs (- (square guess) x)) 0.001))
;; (define (improve guess x)
;;   (average guess (/ x guess)))


;; (define (fixed-point f first-guess)
;;   (define (close-enough? v1 v2)
;;     (display v1)(newline)
;;     (< (abs (- v1 v2))
;;        tolerance))
;;   (define (try guess)
;;     (let ((next (f guess)))
;;       (if (close-enough? guess next)
;;         next
;;         (try next))))
;;   (try first-guess))

(define tolerance 0.00001)
(define (iterative-improve x guess improve)
  (define (good-enough? old-guess guess)
    (< (abs (- old-guess guess)) tolerance))
  (define (iter old-guess guess)
    (if (good-enough? old-guess guess)
      guess
      (iter guess (improve guess))))
  (iter (- guess 20) guess))


(define (fixed-point f first-guess)
  (define (improve guess)
    (f guess))
  (iterative-improve f 1.0 improve))

(define (sqrt x)
  (define (improve guess)
    (average guess (/ x guess)))
  (iterative-improve x 1.0 improve))

(define (logb b x)
  (/ (log x) (log b)))

(define (nrt x n)
  (define (r n) (floor (logb 2 n)))
  (fixed-point-of-transform
    (lambda (y) (/ x (expt y (- n 1)))) (repeated average-damp (r n)) 1.0))

(define (make-rat n d)
  (let ((n (/ (* n d) (abs d)))
        (d (abs d)))
    (let ((g (gcd n d)))
      (cons (/ n g) (/ d g)))))

(define (numer x) (car x))
(define (denom x) (cdr x))

(define (print-rat x)
  (display (numer x))
  (display "/")
  (display (denom x)))

(define (make-point x y)
  (cons x y))

(define (make-segment start end)
  (cons start end))

(define (midpoint segment)
  (make-point
    (average (point-x (car segment)) (point-x (cdr segment)))
    (average (point-y (car segment)) (point-y (cdr segment)))))

; up-left point length width
(define (make-rectangle down-left-point l w)
  (list 0 down-left-point l w))

; diagonal and angle
(define (make-rectangle-2 down-left-point l w)
  (list 1
        down-left-point
        (sqrt (+ (square l)
                 (square w)))
        (atan (/ w l))))

(define (len r)
  (cond ((= (car r) 0) (caddr r))
        ((= (car r) 1) (* (caddr r) (cos (cadddr r))))))

(define (width r)
  (cond ((= (car r) 0) (cadddr r))
        ((= (car r) 1) (* (caddr r) (sin (cadddr r))))))

(define (perimeter r)
  (* 2 (+ (len r) (width r))))

(define (area r)
  (* (len r) (width r)))

(define (point-x p)
  (car p))

(define (point-y p)
  (cdr p))

;; (define (cons x y)
;;   (lambda (m) (m x y)))
;;
;; (define (car z)
;;   (z (lambda (p q) p)))
;;
;; (define (cdr z)
;;   (z (lambda (p q) q)))

(define zero (lambda (f) (lambda (x) x)))

(define (add-1 n)
(lambda (f) (lambda (x) (f ((n f) x)))))

(define (add-interval x y)
  (make-interval (+ (lower-bound x) (lower-bound y))
                 (+ (upper-bound x) (upper-bound y))))

(define (sub-interval x y)
  (make-interval (- (lower-bound x) (upper-bound y))
                 (- (upper-bound x) (lower-bound y))))

(define (mul-interval x y)
  (let ((p1 (* (lower-bound x) (lower-bound y)))
        (p2 (* (lower-bound x) (upper-bound y)))
        (p3 (* (upper-bound x) (lower-bound y)))
        (p4 (* (upper-bound x) (upper-bound y))))
    (make-interval (min p1 p2 p3 p4)
                   (max p1 p2 p3 p4))))

(define (div-interval x y)
  (if (or (= (lower-bound y) 0) (= (upper-bound y) 0))
    (error "Division by 0"))
  (mul-interval
    x
    (make-interval (/ 1.0 (upper-bound y))
                   (/ 1.0 (lower-bound y)))))

(define (make-interval a b) (cons a b))
(define (lower-bound interval) (car interval))
(define (upper-bound interval) (cdr interval))

(define (make-center-width c w)
  (make-interval (- c w) (+ c w)))

(define (center i)
  (/ (+ (lower-bound i) (upper-bound i)) 2))

(define (width i)
  (/ (- (upper-bound i) (lower-bound i)) 2))

(define (percent i)
  (* (/ (width i) (center i)) 100))

(define (make-center-percent center percent)
  (let ((width (* center percent (/ 1 100))))
    (make-center-width center width)))

(define (last-pair l)
  (list (list-ref l (- (length l) 1))))

(define (reverse l)
  (define (iter i res)
    (if (>= i (length l))
      res
      (iter (+ i 1) (cons (list-ref l i) res))))
  (iter 0 (list)))

(define (first-denomination coin-values)
  (car coin-values))

(define (except-first-denomination coin-values)
  (cdr coin-values))

(define (no-more? coin-values)
  (null? coin-values))

(define (cc amount coin-values)
  (cond ((= amount 0) 1)
        ((or (< amount 0) (no-more? coin-values)) 0)
        (else
          (+ (cc amount
                 (except-first-denomination
                   coin-values))
             (cc (- amount
                    (first-denomination
                      coin-values))
                 coin-values)))))

(define (same-parity . args)
  (let ((x (car args)))
    (define (iter i res)
      (if (>= i (length args))
        res
        (let ((y (list-ref args i)))
          (if (= (remainder x 2) (remainder y 2))
            (iter (+ i 1) (append res (list y)))
            (iter (+ i 1) res)))))
    (iter 0 '())))

;; (define (square-list items)
;;   (if (null? items)
;;     nil
;;     (cons (square (car items))
;;           (square-list (cdr items)))))

(define (square-list items)
  (map square items))

(define (for-each f items)
  (unless (null? items)
    (begin (f (car items))
           (for-each f (cdr items)))))


;; (define (square-list items)
;;   (define (iter things answer)
;;     (if (null? things)
;;       answer
;;       (iter (cdr things)
;;             (cons (square (car things))
;;                   answer))))
;;   (iter items nil))
;;
;; (define (square-list items)
;;   (define (iter things answer)
;;     (if (null? things)
;;       answer
;;       (iter (cdr things)
;;             (append answer
;;                   (list (square (car things)))))))
;;   (iter items nil))

(define (deep-reverse l)
  (define (iter i res)
    (if (>= i (length l))
      res
      (let ((e (list-ref l i)))
        (iter
          (+ i 1)
          (cons (if (list? e)
                  (deep-reverse e) e)
                res)))))
  (iter 0 (list)))

(define (count-leaves x)
  (cond ((null? x) 0)
        ((not (pair? x)) 1)
        (else (+ (count-leaves (car x))
                 (count-leaves (cdr x))))))

(define (fringe x)
  (cond ((null? x) x)
        ((not (pair? x)) (list x))
        (else (append
                (fringe (car x))
                (fringe (cdr x))))))

(define (make-mobile left right)
  (cons left right))

(define (make-branch length structure)
  (cons length structure))

(define (left-branch mobile)
  (car mobile))

(define (right-branch mobile)
  (cdr mobile))

(define (branch-length branch)
  (car branch))

(define (branch-structure branch)
  (cdr branch))

(define (is-mobile? s)
  (pair? s))

(define (weight s)
  (if (is-mobile? s) (total-weight s) s))

(define (total-weight x)
  (define (branch-weight b)
    (let ((s (branch-structure b)))
      (weight s)))
  (let ((l (left-branch x))
        (r (right-branch x)))
    (+ (branch-weight l) (branch-weight r))))

(define (balanced? m)
  (if (not (is-mobile? m))
    #t
    (let ((l (left-branch m))
          (r (right-branch m)))
      (and
        (=
          (* (branch-length l) (weight (branch-structure l)))
          (* (branch-length r) (weight (branch-structure r))))
        (balanced? (branch-structure l))
        (balanced? (branch-structure r))))))

(define a
  (make-mobile
    (make-branch 5 10)
    (make-branch 2 (make-mobile
                     (make-branch 5 (make-mobile (make-branch 1 2)
                                                 (make-branch 3 4)))
                     (make-branch 8 9)))))

(define b
  (make-mobile
    (make-branch 5 10)
    (make-branch 2
                 (make-mobile
                   (make-branch 1 2)
                   (make-branch 3 4)))))

(define bb
  (make-mobile
    (make-branch 1 b)
    (make-branch 1 b)))


(define bt
  (make-mobile
    (make-branch 8 2)
    (make-branch 4 4)))

(define bbt
  (make-mobile
    (make-branch 1 bt)
    (make-branch 1 bt)))

(define (scale-tree tree factor)
  (map (lambda (sub-tree)
         (if (pair? sub-tree)
           (scale-tree sub-tree factor)
           (* sub-tree factor)))
       tree))

(define (map-tree f tree)
  (map (lambda (sub-tree)
         (if (pair? sub-tree)
           (map-tree f sub-tree)
           (f sub-tree)))
       tree))

(define (square-tree tree)
  (map-tree square tree))

(define (subsets s)
  (if (null? s)
    (list nil)
    (let ((rest (subsets (cdr s))))
      (append rest (map (lambda (x) (cons (car s) x)) rest)))))

(define (accumulate op initial sequence)
  (if (null? sequence)
    initial
    (op (car sequence)
        (accumulate op initial (cdr sequence)))))

;; (define (map p sequence)
;;   (accumulate (lambda (x y) (cons (p x) y)) nil sequence))
;;
;; (define (append seq1 seq2)
;;    (accumulate cons seq2 seq1))
;;
;; (define (length sequence)
;;   (accumulate (lambda (x y) (+ y 1)) 0 sequence))

(define (horner-eval x coefficient-sequence)
  (accumulate (lambda (this-coeff higher-terms)
                (+ this-coeff (* higher-terms x)))
              0
              coefficient-sequence))

(define (count-leaves x)
  (cond ((null? x) 0)
        ((not (pair? x)) 1)
        (else (+ (count-leaves (car x))
                 (count-leaves (cdr x))))))

(define (count-leaves t)
  (accumulate + 0
    (map (lambda (node)
           (cond ((null? node) 0)
                 ((pair? node) (count-leaves node))
                 (else 1))) t)))

(define (accumulate-n op init seqs)
  (if (null? (car seqs))
    nil
    (cons (accumulate op init (map car seqs))
          (accumulate-n op init (map cdr seqs)))))

;; (define (matrix-*-vector m v)
;;   (map
;;     (lambda (row)
;;       ; TODO: can we use map here?
;;       ; e.g. (map + '(1 2 3) '(3 2 1)) gives (4 4 4)
;;       ; then you can accumulate (4 4 4) with + into 12
;;       (define (iter i res)
;;         (if (>= i (length row))
;;           res
;;           (iter (+ i 1)
;;                 (+ res
;;                    (* (list-ref row i)
;;                       (list-ref v i))))))
;;       (iter 0 0))
;;     m))

;; e.g. if you have
;; (1 2 3)
;; (4 5 6)
;; One row is (1 2 3) to multiply by e.g. (3 4) and (9 0)
;;
;;
;;

; this is the scalar product of the vectors
; eg (1 2 3) * (1 2 3) => (1*1) + (2*2) + (3*3)
(define (dot-product v w)
  (accumulate + 0 (map * v w)))

(define (matrix-*-vector m v)
  (map (lambda (row) (accumulate + 0 (map * row v))) m))

(define (transpose mat)
  (accumulate-n cons nil mat))

(define (matrix-*-matrix m n)
  (define (row-for-col row col)
    (accumulate + 0 (map * row col)))
  (let ((cols (transpose n)))
    (map (lambda (row)
      (map (lambda (col) (row-for-col row col)) cols)) m)))

(define fold-right accumulate)

(define (fold-left op initial sequence)
  (define (iter result rest)
    (if (null? rest)
      result
      (iter (op result (car rest))
            (cdr rest))))
  (iter initial sequence))

(define (reverse sequence)
  (fold-left (lambda (x y) (cons y x)) nil sequence))

(define (reverse sequence)
  (fold-right (lambda (x y) (append y (list x))) nil sequence))

(define (enumerate-interval low high)
  (if (> low high)
    nil
    (cons low (enumerate-interval (+ low 1) high))))

(define (flatmap proc seq)
  (accumulate append nil (map proc seq)))

(define (prime-sum? pair)
  (prime? (+ (car pair) (cadr pair))))

; WARNING: pair here is (1 2), not (1 . 2)
; that is (list 1 2), not (cons 1 2)
(define (make-pair-sum pair)
  (list (car pair) (cadr pair) (+ (car pair) (cadr pair))))

(define (prime-sum-pairs n)
  (map make-pair-sum
       (filter prime-sum? (unique-pairs n))))

(define (permutations s)
  (if (null? s) ; empty set?
    (list nil) ; sequence containing empty set
    (flatmap (lambda (x)
               (map (lambda (p) (cons x p))
                    (permutations (remove x s))))
             s)))

(define (remove item sequence)
  (filter (lambda (x) (not (= x item)))
          sequence))

;; this is adapted from the book code for (prime-sum-pairs n)
;; e.g. for 3 -> (2 1) (3 1) (3 2)
;;
;; (define (unique-pairs n)
;;   (flatmap
;;     (lambda (i)
;;       (map (lambda (j) (list i j))
;;            (enumerate-interval 1 (- i 1))))
;;     (enumerate-interval 1 n)))

;; my way: this produces ordered pairs, e.g. for 3 -> (1 2) (1 3) (2 3)
;;
(define (unique-pairs n)
  (flatmap (lambda (i)
         (map (lambda (j) (list i j))
              (enumerate-interval (+ i 1) n)))
       (enumerate-interval 1 n)))

;; same for triples, just add a nested lambda
(define (unique-triples n)
  (flatmap
    (lambda (i)
      (flatmap
        (lambda (j)
          (map (lambda (k) (list i j k))
               (enumerate-interval (+ j 1) n)))
        (enumerate-interval (+ i 1) n)))
    (enumerate-interval 1 n)))

;; general implementation for k-tuples
;; (define (unique-tuples n k)
;;      (define (iter m k)
;;          (if (= k 0)
;;              (list nil)
;;              (flatmap (lambda (j)
;;                          (map (lambda (tuple) (cons j tuple))
;;                              (iter (+ j 1) (- k 1))))
;;                      (enumerate-interval m n))))
;;      (iter 1 k))

; define pairs and triples in terms of k-tuples
;; (define (unique-pairs n) (unique-tuples n 2))
;; (define (unique-triples n) (unique-tuples n 3))

(define (unique-triples-sum n s)
  (filter
    (lambda (e) (= (accumulate + 0 e) s))
    (unique-triples n)))


(define (make-position row col)
  (cons row col))

(define (position-row pos)
  (car pos))

(define (position-col pos)
  (cdr pos))

(define empty-board '())

(define (adjoin-position row col positions)
  (append positions (list (make-position row col))))

(define (safe? col positions)
   (let ((kth-queen (list-ref positions (- col 1)))
         (other-queens (filter (lambda (q)
                                 (not (= col (position-col q))))
                               positions)))
   (define (attacks? q1 q2)
     (or (= (position-row q1) (position-row q2))
         (= (abs (- (position-row q1) (position-row q2)))
            (abs (- (position-col q1) (position-col q2))))))

   (define (iter q board)
     (or (null? board)
         (and (not (attacks? q (car board)))
              (iter q (cdr board)))))
   (iter kth-queen other-queens)))


(define (queens board-size)
  (define (queen-cols k)
    (if (= k 0)
      (list empty-board)
      (filter
        (lambda (positions) (safe? k positions))
        (flatmap
          (lambda (rest-of-queens)
            (map (lambda (new-row)
                   (adjoin-position
                     new-row k rest-of-queens))
                 (enumerate-interval 1 board-size)))
          (queen-cols (- k 1))))))
  (queen-cols board-size))

(define true #t)
(define false #f)

(define (equal? l1 l2)
  (cond ((and (null? l1) (null? l2)) true)
        ((or
           (not (eq? (null? l1) (null? l2)))
           (not (eq? (symbol? l1) (symbol? l2))) false))
        ((and (not (list? l1)) (not (list? l2))) (eq? l1 l2))
        ((eq? (car l1) (car l2)) (equal? (cdr l1) (cdr l2)))
        (else false)))

(define (deriv exp var)
  (cond ((number? exp) 0)
        ((variable? exp) (if (same-variable? exp var) 1 0))
        ((sum? exp) (make-sum (deriv (addend exp) var)
                              (deriv (augend exp) var)))
        ((product? exp)
         (make-sum
           (make-product (multiplier exp)
                         (deriv (multiplicand exp) var))
           (make-product (deriv (multiplier exp) var)
                         (multiplicand exp))))
        ((exponentiation? exp)
         (make-product
           (make-product (exponent exp)
                         (make-exponentiation
                           (base exp)
                           (- (exponent exp) 1)))
           (deriv (base exp) var)))
        (else
          (error "unknown expression type: DERIV" exp))))

(define (variable? x) (symbol? x))

(define (same-variable? v1 v2)
  (and (variable? v1) (variable? v2) (eq? v1 v2)))

(define (=number? exp num) (and (number? exp) (= exp num)))

(define (make-sum a1 a2)
  (cond ((=number? a1 0) a2)
        ((=number? a2 0) a1)
        ((and (number? a1) (number? a2))
         (+ a1 a2))
        (else (list a1 '+ a2))))

(define (make-product m1 m2)
  (cond ((or (=number? m1 0) (=number? m2 0)) 0)
        ((=number? m1 1) m2)
        ((=number? m2 1) m1)
        ((and (number? m1) (number? m2)) (* m1 m2))
        (else (list m1 '* m2))))

(define (make-exponentiation b e)
  (cond ((=number? b 0) 0)
        ((=number? e 0) 1)
        ((=number? e 1) b)
        ((and (number? b) (number? e)) (expt b e))
        (else (list b '^ e))))

(define (sum? x) (and (pair? x) (eq? (cadr x) '+)))
(define (product? x) (and (pair? x) (eq? (cadr x) '*)))
(define (exponentiation? x) (and (pair? x) (eq? (cadr x) '^)))

(define (addend s) (car s))
(define (augend s) (accumulate make-sum 0 (cddr s)))

(define (multiplier p) (car p))
(define (multiplicand p) (accumulate make-product 1 (cddr p)))

(define (base s) (car s))
(define (exponent s) (caddr s))

;; set as unordered lists - no duplicates
(define (element-of-set? x set)
  (cond ((null? set) false)
        ((equal? x (car set)) true)
        (else (element-of-set? x (cdr set)))))

(define (adjoin-set x set)
  (if (element-of-set? x set)
    set
    (cons x set)))

(define (intersection-set set1 set2)
  (cond ((or (null? set1) (null? set2)) '())
        ((element-of-set? (car set1) set2)
         (cons (car set1) (intersection-set (cdr set1) set2)))
        (else (intersection-set (cdr set1) set2))))

(define (union-set set1 set2)
  (cond ((null? set1) set2)
        ((null? set2) set1)
        (else (adjoin-set (car set1)
                          (union-set (cdr set1) set2)))))

;; sets as unordered lists - duplicates
;; (define (adjoin-set x set)
;;   (cons x set))
;;
;; (define (union-set set1 set2)
;;   (append set1 set2))

;; sets as ordered list
(define (element-of-set? x set)
  (cond ((null? set) false)
        ((= x (car set)) true)
        ((< x (car set)) false)
        (else (element-of-set? x (cdr set)))))

(define (intersection-set set1 set2)
  (if (or (null? set1) (null? set2))
    '()
    (let ((x1 (car set1)) (x2 (car set2)))
      (cond ((= x1 x2)
             (cons x1 (intersection-set (cdr set1)
                                        (cdr set2))))
            ((< x1 x2)
             (intersection-set (cdr set1) set2))
            ((< x2 x1)
             (intersection-set set1 (cdr set2)))))))

(define (adjoin-set x set)
  (cond ((element-of-set? x set) set)
        ((null? set) (list x))
        ((< x (car set)) (cons x set))
        (else (cons (car set) (adjoin-set x (cdr set))))))

(define (union-set set1 set2)
  (cond ((null? set1) set2)
        ((null? set2) set1)
        (else
          (let ((x1 (car set1))
                (x2 (car set2)))
            (cons (min x1 x2)
                  (union-set (if (> x1 x2)
                                 set1
                                 (cdr set1))
                               (if (> x2 x1)
                                 set2
                                 (cdr set2))))))))

; tree representation of sets
(define (entry tree) (car tree))
(define (left-branch tree) (cadr tree))
(define (right-branch tree) (caddr tree))
(define (make-tree entry left right)
  (list entry left right))

(define (adjoin-set x set)
  (cond ((null? set) (make-tree x '() '()))
        ((= x (entry set)) set)
        ((< x (entry set))
         (make-tree (entry set)
                    (adjoin-set x (left-branch set))
                    (right-branch set)))
        ((> x (entry set))
         (make-tree (entry set) (left-branch set)
                    (adjoin-set x (right-branch set))))))

(define (tree->list-1 tree)
  (if (null? tree)
    '()
    (append (tree->list-1 (left-branch tree))
            (cons (entry tree)
                  (tree->list-1
                    (right-branch tree))))))

(define (tree->list-2 tree)
  (define (copy-to-list tree result-list)
    (if (null? tree)
      result-list
      (copy-to-list (left-branch tree)
                    (cons (entry tree)
                          (copy-to-list
                            (right-branch tree)
                            result-list)))))
  (copy-to-list tree '()))

