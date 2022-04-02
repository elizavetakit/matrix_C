SRCS = s21_matrix.c
	
OBJS = ${SRCS:.c=.o}
CC = gcc
CFLAGS = -Wall -Werror -Wextra -std=c11 -D_GNU_SOURCE -lm

all: clean s21_matrix.a test

s21_matrix.a: ${SRCS}
	${CC} ${CFLAGS} -c ${SRCS}
	@ar rcs s21_matrix.a ${OBJS}
	@ranlib s21_matrix.a

test: s21_matrix.a unit_test

unit_test:
	${CC} ${CFLAGS} -c s21_matrix_tests.c
	${CC} s21_matrix_tests.o s21_matrix.a -lcheck -o test
	./test

gcov_report: s21_matrix.a
	@gcc --coverage s21_matrix_tests.c -lcheck ${SRCS} -o unit-test
	@./unit-test
	@lcov -t "test" -o test.info -c -d .
	@genhtml -o report test.info
	
clean:
	@/bin/rm -rf *.o *.a test unit-test *.gcno *gcda report *.info main

checks: leaks style

style:
	python3 ../materials/linters/cpplint.py --extensions=c *.c *.h

cppcheck:
	cppcheck --enable=all --suppress=missingIncludeSystem . *.c

leaks: s21_matrix_tests.c
	CK_FORK=no leaks --atExit -- ./test
	