echo "task2_test_camera ..."
../../build/examples/task1/task2_test_camera_model

echo "task3_test_matching ..."
../../build/examples/task1/task3_test_feature_matching ../data/kxm1.jpg ../data/kxm2.jpg

echo "task4_test_fundamental_8_point ..."
../../build/examples/task1/task4_test_fundamental_8_point

echo "task5_test_fundamental_ransac ..."
../../build/examples/task1/task5_test_fundamental_ransac ./correspondences.txt

echo "task6_test_pose_from_fundamental ..."
../../build/examples/task1/task6_test_pose_from_fundamental

