NAME		=	ConsoleResonator

CLASSES		=	Resonator Linear_1D Linear_2D NonLinear_1D
FILES		=	$(NAME) $(CLASSES)
SRCS		=	$(patsubst %,%.cpp,$(FILES))
OBJS		=	$(patsubst %,%.o,$(FILES))
HEADERS		=	$(patsubst %,%.h,$(CLASSES))

CPP			=	g++
FLAGS		=	-Wall -Wextra -Werror -fopenmp -O3
RM			=	/bin/rm -f

.PHONY: all clean fclean re

all:		$(NAME)

$(NAME):	$(OBJS) $(HEADERS)
			$(CPP) $(OBJS) $(FLAGS) -o $(NAME)

clean:
			$(RM) $(OBJS)

fclean:		clean
			$(RM) $(NAME)

re:			fclean	all
