#pragma once
#include "solve.h"
#include "test.h"
#include <Windows.h>

int main() {
    SetConsoleOutputCP(1251);
    srand(GetTickCount64());
    

    /*const std::vector<double> w = { 1, 2, 3};
    const std::vector<double> lambdas = {-14, 21, 17};
    std::vector<std::vector<double>> H;*/

    //in.A = CreateA(w, lambdas, H);
    
    //testsRun();

	std::cout << "1. Чтение из файла\n2. Запустить тестирование\n";
	int ch;
	std::cin >> ch;
	switch (ch)
	{
	case 1:
	{
		InParams in;
		OutParams out;
		std::string filename = "input.txt";
		ReadInputFromFile(filename, in);
		Solve(in, out);
		print(out);
		break;
	}
	case 2:
	{
		testsRun();
		break;
	}
	default:
		std::cout << "Ошибка ввода\n";
	}

    std::cin.ignore().get();
    return 0;
}
