#ifndef READ_FILE_TO_ARRAY																//read_file_to_array.h
#define READ_FILE_TO_ARRAY

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

namespace READ_FILE_TO_ARRAY
{
	using namespace std;
//возвращает файл с записанными градиентами
	int read_files_for_elelastic (	dealii::Vector<dbl> &x_coordinate_of_vertex,				//ссылка на векторы с координатами
						dealii::Vector<dbl> &y_coordinate_of_vertex,							//ссылка на векторы с координатами
						dealii::Vector<dbl> &U_x, dealii::Vector<dbl> &U_y,						//ссылки на векторы значений
						dealii::Vector<dbl> &U_x_gradX, dealii::Vector<dbl> &U_x_gradY,			//ссылка на векторы с градиентами
						dealii::Vector<dbl> &U_y_gradX, dealii::Vector<dbl> &U_y_gradY,			//ссылка на векторы с градиентами
						int Colons)																//Количество столбцов в файле
	{
		string line;																//сохраняет строку из файла
		string subline;																//сохраняет строку из файла
		double number;
		ifstream file_in ("out/GRAD.gpl");											//поток файла
		int StrNum = 0;																//инициализация индекса для массива
		int ColNum = 0;																//инициализация столбца
		if (file_in.is_open())
		{
			while ( getline (file_in,line) )										//пока есть строки в файле
			{
				for (int i=0; i<line.size(); ++i)									//Пробежать посимвольно в строке
				{
					if(line[i] != ' ' && line[i] != '\t' && i < line.size()-1 )		//Проверить является ли символ числом и добавить в подстроку
					{
						subline += line[i];
					}
					else															//Если подстрока выделена, то вывести её и обнулить
					{
						if (subline != "")
						{
							number = atof(subline.c_str());
//							cout << number << "\t\t";								//Вывод числа
							switch(ColNum)
							{
								case 0:
									{												//Будут содержаться значения номера строки, начиная с нуля
										StrNum = number;
									}
								case 1:
									{												//Будут содержаться значения x_coordinate_of_vertex
										!(StrNum % 2) ?	  x_coordinate_of_vertex(StrNum / 2) = number						//чётное
														: x_coordinate_of_vertex(StrNum / 2) = number;						//нечётное
									}
								case 2:
									{												//Будут содержаться значения y_coordinate_of_vertex
										!(StrNum % 2) ?	  y_coordinate_of_vertex(StrNum / 2) = number						//чётное
														: y_coordinate_of_vertex(StrNum / 2) = number;						//нечётное
									}
								case 3:
									{												//Будут содержаться значения U_x и U_y
										!(StrNum % 2) ?	  U_x(StrNum / 2) = number						//чётное
														: U_y(StrNum / 2) = number;						//нечётное
									}
								case 4:
									{												//Будут содержаться производные по x
										!(StrNum % 2) ?	  U_x_gradX(StrNum / 2) = number						//чётное
														: U_y_gradX(StrNum / 2) = number;						//нечётное
									}
								case 5:
									{												//Будут содержаться производные по y
										!(StrNum % 2) ?	  U_x_gradY(StrNum / 2) = number						//чётное
														: U_y_gradY(StrNum / 2) = number;						//нечётное
//										cout << "massivY(" << NumberOfStringInFile << ") = " << number << "\t\t";					//Вывод полученного числа
									}
							}
							if(ColNum != Colons-1) {ColNum++;}		//увеличить индекс. Смотрим с какой колонки мы считали число
							else {ColNum = 0;}										//увеличить индекс. Смотрим с какой колонки мы считали число
						}
						subline = "";
					}
				}
//				NumberOfStringInFile++;												//увеличить индекс. К номеру строки прибавим 1
//				cout << "\n";														//Перевод на новую сроку после окончания обработки одной строки из файла
			}
			file_in.close();
		}
		else 
		{
			cout << "READ_FILE_TO_ARRAY::read_files:   Unable to open file";
			file_in.close();
		}
		return 0;
	}


//возвращает файл с записанными градиентами
int read_files_for_temperature (	dealii::Vector<dbl> &x_coordinate_of_vertex,		//ссылка на векторы с координатами
						dealii::Vector<dbl> &y_coordinate_of_vertex,					//ссылка на векторы с координатами
						dealii::Vector<dbl> &U_z,										//ссылка на вектор значений
						dealii::Vector<dbl> &U_z_gradX,									//ссылка на вектор с градиентом
						dealii::Vector<dbl> &U_z_gradY,									//ссылка на вектор с градиентом
						int Colons)														//Количество столбцов в файле
	{
		string line;																//сохраняет строку из файла
		string subline;																//сохраняет строку из файла
		double number;
		ifstream file_in ("out/GRAD.gpl");											//поток файла
		int StrNum = 0;																//инициализация индекса для массива
		int ColNum = 0;																//инициализация столбца
		if (file_in.is_open())
		{
			while ( getline (file_in,line) )										//пока есть строки в файле
			{
				for (int i=0; i<line.size(); ++i)									//Пробежать посимвольно в строке
				{
					if(line[i] != ' ' && line[i] != '\t' && i < line.size()-1 )		//Проверить является ли символ числом и добавить в подстроку
					{
						subline += line[i];
					}
					else															//Если подстрока выделена, то вывести её и обнулить
					{
						if (subline != "")
						{
							number = atof(subline.c_str());
//							cout << number << "\t\t";								//Вывод числа
							switch(ColNum)
							{
								case 0:
									{												//Будут содержаться значения номера строки, начиная с нуля
										StrNum = number;
									}
								case 1:
									{												//Будут содержаться значения x_coordinate_of_vertex
										x_coordinate_of_vertex(StrNum) = number;
									}
								case 2:
									{												//Будут содержаться значения y_coordinate_of_vertex
										y_coordinate_of_vertex(StrNum) = number;
									}
								case 3:
									{												//Будут содержаться значения U_z
										U_z(StrNum) = number;
									}
								case 4:
									{												//Будут содержаться производные по x
										U_z_gradX(StrNum) = number;
									}
								case 5:
									{												//Будут содержаться производные по y
										U_z_gradY(StrNum) = number;
									}
							}
							if(ColNum != Colons-1) {ColNum++;}		//увеличить индекс. Смотрим с какой колонки мы считали число
							else {ColNum = 0;}										//увеличить индекс. Смотрим с какой колонки мы считали число
						}
						subline = "";
					}
				}
//				NumberOfStringInFile++;												//увеличить индекс. К номеру строки прибавим 1
//				cout << "\n";														//Перевод на новую сроку после окончания обработки одной строки из файла
			}
			file_in.close();
		}
		else cout << "READ_FILE_TO_ARRAY::read_files_for_temperature:   Unable to open file";
		file_in.close();
		return 0;
	}



//записывает в файл 
void read_vector(const char *FileName,						//Строка с названием файла и его расположением
					dealii::Vector<dbl> &AnyVec,			//содержит ссылку на вектор, куда считываются данные
					int VecCol,								//Номер нужного для чтения столбца, нумерация с нуля
					int Colons)								//количество строк
	{
		string line;																//сохраняет строку из файла
		string subline;																//сохраняет строку из файла
		double number;
//		ifstream file_in ("out/GRAD.gpl");											//поток файла
//		ifstream file_in ("out/FU_x_grad.gpl");										//поток файла
//		ifstream file_in ("out/FU_x.gpl");											//поток файла
		ifstream file_in (FileName);
		int StrNum = 0;																//инициализация индекса для массива
		int ColNum = 0;																//инициализация столбца
int nn = 0;
		if (file_in.is_open())
		{
//								std::cout << "\n\nColNum = " << ColNum << "\n";
//								std::cout << "file_in = " << file_in << "\n";
//								std::cout << "line    = " << line << "\n";
//								std::cout << "getline (file_in,line) = " << getline (file_in,line) << "\n\n";
//			while ( !file_in.eof() )										//пока есть строки в файле
			for (nn = 0; nn < AnyVec.size(); nn++)
			{
				getline (file_in,line);
				for (int i=0; i<line.size(); ++i)									//Пробежать посимвольно в строке
				{
					if(line[i] != ' ' && line[i] != '\t' && i < line.size()-1 )		//Проверить является ли символ числом и добавить в подстроку
					{
						subline += line[i];
					}
					else															//Если подстрока выделена, то вывести её и обнулить
					{
						if (subline != "")
						{
							number = atof(subline.c_str());
//							cout << number << "\t222\t\n";
							if(ColNum == 0)
							{
								StrNum = number;
							}
							if(ColNum == VecCol)
							{
								AnyVec[StrNum] = number;
//								std::cout << "AnyVec[" << StrNum << "] = " << AnyVec[StrNum] << "\t\tnumber = " << number << "\n";
//								std::cout << "qqqqqq\n";
//								std::cout << "ColNum = " << ColNum << "\n";
							}
							if(ColNum != Colons-1) {ColNum++;}						//увеличить индекс. Смотрим с какой колонки мы считали число
							else {ColNum = 0;}										//увеличить индекс. Смотрим с какой колонки мы считали число
						}
						subline = "";
					}
				}
//				NumberOfStringInFile++;												//увеличить индекс. К номеру строки прибавим 1
//				cout << "\n";														//Перевод на новую сроку после окончания обработки одной строки из файла
			}
			file_in.close();
		}
		else cout << "READ_FILE_TO_ARRAY::read_files_for_temperature:   Unable to open file";
		file_in.close();


//		return 0;
	}






}
#endif
