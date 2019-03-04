git checkout develop && 
git add SbgatCore/include && 
git add SbgatCore/source && 
git add Examples/ && 
git commit -m "shipping to fortuna" &&
git push origin develop &&
git checkout fortuna && 
git checkout develop SbgatCore/include && 
git checkout develop SbgatCore/source && 
git checkout develop 'Examples/**/*.py' &&
git checkout develop 'Examples/**/*.cpp' &&
git commit -m "shipping to fortuna" &&
git push origin fortuna &&
git checkout develop
