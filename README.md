---
layout: post
title: Poisson Denklemi Yardımıyla Görüntü Düzenleme
slug: poisson-image-editing
author: Bahri ABACI
categories:
- Görüntü İşleme Uygulamaları
- Lineer Cebir
- Nümerik Yöntemler
thumbnail: resources//thumbnail.png 
---
Görüntü düzenleme (image editing), imge üzerinde yer alan renklerin evrensel veya yerel ilişkiler ile değiştirilmesi işlemidir. Önceki yazılarımızda bahsettiğimiz Renk Dönüşümleri, Instagram Filtreleme Yöntemleri gibi yöntemler imge üzerinde evrensel ilişkiler ile dönüşüm yaparken, bu yazımızın konusu olan Poisson eşitliği yardımıyla görüntü düzenleme algoritması yerel ilişkilere dayalı bir dönüşüm işlemidir. Patrick Pérez ve arkadaşları tarafından 2003 yılında önerilen bu yöntem basit bir anlatımla verilen bir imgenin (kaynak:source) gradyanını elimizde bulunan bir başka imgeye (hedef:target) kopyalama işlemidir. Bu kopyalama işlemi sonucunda kaynak imgenin dokusu (gradyanı) hedef imge içerisine yerleştirilmiş olur.

<!--more-->
$\mathbf{I}$ çıktı imgesini, $\mathbf{T}$ kopyalamanın yapılacağı arka plan imgesini, $\mathbf{S}$ kopyalanacak kaynak imgeyi ve $\Omega$ kopyalanacak bölgeye ait sınır maskesini göstermek üzere; önerilen yöntem aşağıdaki hata fonksiyonunu en küçüklemeye çalışmaktadır.

$$
\begin{aligned}
\min_{I(x,y) \in \Omega} \sum_{\Omega} \lVert\nabla I(x,y)-\nabla S(x,y)\rVert \\
\text{subject to} \quad I(x,y)=T(x,y) \text{ on } \partial \Omega
\end{aligned}
\label{cost} \tag{1}
$$

Denklem \ref{cost} dan da görüldüğü üzere yöntem hata değeri olarak $\Omega$ bölgesi içerisinde çıktının imgesinin gradyanını kaynak imgenin gradyanına olabildiğince yakın olmaya zorlamaktadır. Tek başına bu maliyet fonksiyonu $I(x,y) = S(x,y),\; \forall x,y \in \Omega$ olmayı zorlasa da denklemde yer alan $I(x,y)=T(x,y), \; \forall x,y \in \partial \Omega$ kısıtı kaynak imgenin izsiz (seamless) bir şekilde hedef imgeye yerleştirilmesini sağlamaktadır. Yöntemin anlaşılması için matematiksel detaylarına devam etmeden önce şu ana kadar verdiğimiz matematiksel ifadeleri bir örnek üzerinde inceleyelim.

| $\mathbf{T}$: Hedef İmge |  $\mathbf{S}$: Kaynak İmge | $\Omega$: Kaynak Maskesi | $\partial \Omega$: Kaynak Maskesi Sınırları | $\mathbf{I}$: Birleştirilmiş İmge
:-------:|:----:|:----:|:---:|:---:|
![Poisson Image Editing][bear_background] | ![Poisson Image Editing][bear_foreground] | ![Poisson Image Editing][bear_foreground_mask] | ![Poisson Image Editing][bear_foreground_mask_explained] | ![Poisson Image Editing][bear_blended_image] |

Yukarıda verilen imgelerde ilk imge formüllerde $\mathbf{T}$ ile ifade edilen ve kopyalamanın yapılacağı arka plan imgesini göstermektedir. $\mathbf{S}$ kaynak imgesi $\mathbf{T}$ imgesi üzerine yerleştirilmek istenen ön plan imgesini göstermektedir. $\Omega$ ile verilen imgede yer alan beyaz bölgeler ise $\mathbf{S}$ imgesinin hangi piksellerinin $\mathbf{T}$ imgesine yerleştirileceğini gösteren maskedir. $\partial \Omega$ ile verilen imge $\Omega$ maskensinin sınır bölgesini göstermektedir. Son olarak $\mathbf{I}$; $\mathbf{T}$ ve $\mathbf{S}$ imgelerinin $\Omega$ maskesi kullanılarak izsiz bir şekilde birleştirilmesi sonucu elde edilen imgeyi göstermektedir.

Denklem \ref{cost} ile verilen maliyet fonksiyonun en küçüklenmesi için ifadeninin gradyanı bulunup sıfıra eşitlenirse; $\nabla \cdot \nabla \mathbf{I} = \nabla \cdot \nabla \mathbf{S}$ olması gerektiği bulunur. Burada $\Delta F = \nabla \cdot \nabla F$, Laplace operatörü olarak bilinmektedir ve iki boyutlu imgeler için 

$$\Delta = \begin{bmatrix}\phantom{+}0 & \phantom{+}1 & \phantom{+}0\\\phantom{+}1 &-4 & \phantom{+}1\\\phantom{+}0 & \phantom{+}1 & \phantom{+}0\end{bmatrix}$$ 

çekirdeği ile evrişim işlemi sonucunda hesaplanır. Burada çözüm hesaplanırken $\partial \Omega$ sınır bölgesinde ise (beyaz bölge) $I(x,y)=T(x,y)$ eşitliğinin sağlanması gerektiği unutulmamalıdır.

Ayrık imgeler için evrişim teorisi kullanılarak $\Delta \mathbf{I} = \Delta \mathbf{S}$ çözümü açık bir şekilde yazılırsa;

$$
\begin{aligned}
4I(x,y)-I(x-1,y)-I(x,y-1)-I(x+1,y)-I(x,y+1)=& \Delta S(x,y) & \text{in} \quad \Omega\\
I(x,y)=&T(x,y) & \text{on} \quad \partial \Omega
\end{aligned}
\label{solution}
\tag{2}
$$

elde edilir. Bulunan ifadede kaynak imge girdi olarak verildiğinden $\Delta S(x,y)$ nin hesaplanmasında herhangi bir sıkıntı bulunmamaktadır. $\mathbf{I}$ imgesi ise yukarıda verilen denklem takımı Poisson matrisi ile aşağıdaki şekilde ifade edilerek bulunulabilir.

$$
\underbrace{
\begin{bmatrix}
&&&&\dots\\
&&&&\dots\\
&&&&\dots\\
\hline
&&&&\dots\\
0 \dots 0 & -1 & 0 \dots 0 & -1 & 4 & -1 & 0 \dots 0 & -1 & 0 \dots 0\\
&&&&\dots\\
\hline
&&&&\dots \\
&&&&\dots \\
&&&&\dots \\
\end{bmatrix}
}_{\mathbf{P}}

\underbrace{
\begin{bmatrix}
\dots \\ I(x,y-1) \\ \dots\\
\hline
I(x-1,y) \\ I(x,y) \\ I(x+1,y)\\ 
\hline
\dots \\I(x,y+1) \\ \dots
\end{bmatrix}
}_{\mathbf{I(\Omega)}}

=

\underbrace{
\begin{bmatrix}
\dots \\ \dots \\ \dots \\
\hline
\dots \\ \Delta S(x,y) \\ \dots \\ 
\hline
\dots \\ \dots \\ \dots
\end{bmatrix}
}_{\Delta \mathbf{S(\Omega)}}
$$

Yukarıda yer alan $\Delta \mathbf{S(\Omega)} = \mathbf{P} \mathbf{I(\Omega)}$ matris ifadesi $\forall x,y,x-1,x+1,y-1,y+1 \; \in \Omega$ olması durumunda geçerlidir. Sınır bölgeleri için $I(x,y) = T(x,y)$ kısıtı bulunduğundan bu noktalar için ifade küçük farklılıklar gösterecektir. Aşağıda örnek maskeden kesilen küçük bir imge parçası için $\forall x,y,x-1,x+1,y-1,y+1 \; \in \Omega$ koşulunun sağlanmaması durumunda denklemin nasıl yazılacağı incelenecektir.

![Poisson Image Editing Border Conditions][poisson_image_editing_border]

Verilen imgede gri bölgeler $\Omega$ imgesindeki siyah (hedef imgeden seçilecek) pikselleri, beyaz (kaynak imgeden seçilecek) pikselleri göstermektedir. Örnek olarak imgede $x,y$ ile işaretlenen piksel için Denklem \ref{solution} ile verilen eşitliği yazarsak;

$$
\begin{aligned}
4I(x,y)-I(x-1,y)-I(x,y-1)-I(x+1,y)-I(x,y+1)=& \Delta S(x,y) & \text{in} \quad \Omega\\
I(x,y)=&T(x,y) & \text{on} \quad \partial \Omega
\end{aligned}
$$

elde edilir. Burada $(x-1,y)$ ve $(x,y-1)$ $\partial \Omega$ sınırında bulunduğundan kısıt gereği $I(x-1,y) = T(x-1,y)$ ve $I(x,y-1) = T(x,y-1)$ olmalıdır. Bu bilgi yukarıdaki eşitlikte yerine yazılırsa;

$$
\begin{aligned}
4I(x,y)-T(x-1,y)-T(x,y-1)-I(x+1,y)-I(x,y+1)=& \Delta S(x,y) \\
4I(x,y)-I(x+1,y)-I(x,y+1)=& \Delta S(x,y) + T(x-1,y) + T(x,y-1)
\end{aligned}
$$

elde edilir. Bu durumda yukarıda tanımlanan matris çarpımı ilgili satır için aşağıdaki şekilde yazılır.

$$
\begin{bmatrix}
&&&&\dots\\
&&&&\dots\\
&&&&\dots\\
\hline
&&&&\dots\\
0 \dots 0 & 0 & 0 \dots 0 & 0 & 4 & -1 & 0 \dots 0 & -1 & 0 \dots 0\\
&&&&\dots\\
\hline
&&&&\dots \\
&&&&\dots \\
&&&&\dots \\
\end{bmatrix}

\begin{bmatrix}
\dots \\ \dots \\ \dots\\
\hline
\dots \\ I(x,y) \\ I(x+1,y)\\ 
\hline
\dots \\I(x,y+1) \\ \dots
\end{bmatrix}

=
\begin{bmatrix}
\dots \\ \dots \\ \dots \\
\hline
\dots \\ \Delta S(x,y) + T(x-1,y) + T(x,y-1)\\ \dots \\ 
\hline
\dots \\ \dots \\ \dots
\end{bmatrix}
$$

$\Omega$ maskesi içerisinde yer alan tüm pikseller için yukarıda verilene benzer şekilde matris eşitlikleri yazıldığında, aranılan $\mathbf{I}$ imgesi $\mathbf{I(\Omega)} = \mathbf{P}^{-1} \left( \mathbf{\Delta S(\Omega)} + \mathbf{T(\partial \Omega)}\right)$ işlemi ile bulunabilir.

Önerilen yöntem ve çözümü oldukça kolay görünse de $\mathbf{P}$ matrisinin boyutu maskede yer alan beyaz piksellerin sayısı ile belirlendiğinden matris tersi bulma işlemi oldukça uzun süreler almaktadır. $\mathbf{P}$ matrisi dikkatli incelendiğinde her satırında en fazla $5$ elamanın sıfırdan farklı olduğu görülür. $\mathbf{P}$ matrisinin bu ayrık yapısı sayesinde matris tersi bulma işlemi Gauss-Seidel veya SOR (Successive Over-Relaxation) gibi iteratif yöntemlerle hızlı bir şekilde hesaplanabilmektedir. Bu yazıda ayrık matrisin tersinin bulunması için Successive Over-Relaxation yöntemi kullanılmıştır.

### Successive Over-Relaxation

SOR yöntemi $A\mathbf {x}=\mathbf {b}$ şeklinde verilen doğrusal denklem takımlarının çözümü için kullanılan bir yöntemdir. Yöntemin çalışma mantığının anlaşılması için bir $A$ kare matrisinin diyagonal $D$, alt üçgen $L$ ve üst üçgen $U$ matrislerinin toplamı şeklinde yazdığımızı düşünelim.

$$
\underbrace{
\begin{bmatrix}
a_{11} & \dots & a_{1n} \\
\vdots & \ddots & \vdots \\ 
a_{n1} & \dots & a_{nn}
\end{bmatrix}
}_{A}
=
\underbrace{
\begin{bmatrix}
a_{11} & \dots & 0 \\
\vdots & \ddots & \vdots \\ 
0 & \dots & a_{nn}
\end{bmatrix}
}_{D}
+
\underbrace{
\begin{bmatrix}
0 & \dots & 0 \\
\vdots & \ddots & \vdots \\ 
a_{n1} & \dots & 0
\end{bmatrix}
}_{L}
+
\underbrace{
\begin{bmatrix}
0 & \dots & a_{1n} \\
\vdots & \ddots & \vdots \\ 
0 & \dots & 0
\end{bmatrix}
}_{U}
$$

$A$ matrisi bu şekilde yazıldıktan sonra $\omega > 1$ esnetme katsayısı kullanılarak $A\mathbf {x}=\mathbf {b}$ ifadesi

$$(D+\omega L)\mathbf {x} = \omega \mathbf {b} - \left[\omega U+(\omega -1)D \right]\mathbf {x}$$ 

şeklinde düzenlenebilir. Burada $(D+\omega L)$ ifadesinin sol tarafa atılarak $\mathbf {x}$ in yalnız bırakılması sonucunda iteratif bir şekilde çalışan SOR yöntemi elde edilmiş olur.

$$\mathbf {x}^{k+1} = (D+\omega L)^{-1} \left( \omega \mathbf {b} - \left[\omega U+(\omega -1)D \right]\mathbf {x}^{k} \right)$$

Elde edilen denklemde yer alan $D+\omega L$ ifadesi alt üçgen bir matris olduğundan tersinin hesaplanması işlemi ileri koyma tekniği ile kolaylıkla yapılabilir. SOR yöntemi bu tekniği de içerecek şekilde güncellenirse çözüm;

$$
x_{i}^{k+1}=(1-\omega )x_{i}^{k}+{\frac {\omega }{a_{ii}}} \left(b_{i}-\sum_{j < i} a_{ij} x_{j}^{k+1} - \sum_{j > i}a_{ij}x_{j}^{k}\right), \; i=1,2,\ldots,n.
$$

eşitliği ile elde edilir. Burada her ne kadar $i=1,2,\ldots,n$ şeklinde yazılsa da $\mathbf{P}$ matrisi ayrık yapıda olduğundan sadece sıfırdan farklı elemanların döngüde kullanılması yeterlidir. $\mathbf{P}$ matrisinin her satırında da en fazla $5$ sıfırdan faklı sayı yer aldığından her bir iterasyon matris boyutundan bağımsız şekilde hesaplanabilmektedir. 

Yazıda anlatılan yöntemin gerçeklenmesi için IMLAB kütüphanesi kullanılmıştır. Kütüphanede ayrık matris tipi bulunmadığından kod içerisinde `struct poisson_element` tipinde bir veri tipi tanımlanmıştır. Bu veri tipi $x,y$ koordinatı için $\mathbf{P}$ matrisinde yer alması gereken katsayıları tutmaktadır. 

$\mathbf{I(\Omega)} = \mathbf{P}^{-1} \left( \mathbf{\Delta S(\Omega)} + \mathbf{T(\partial \Omega)}\right)$ işleminin hesaplanması için SOR başlığı altında anlatılan yöntem `matrix_t* poisson_solver(vector_t *P, matrix_t *S)` fonksiyonunda gerçeklenmiştir. Çalışmada iterasyon sayısı $K=4000$ ve gevşeme katsayısı $\omega=1.8$ olarak kullanılmıştır.

Aşağıda yöntemin farklı imgeler üzerinde çalıştırılması ile elde edilen sonuçlar paylaşılmıştır. Verilen örneklerden de görüldüğü üzere yöntem kaynak imgenin doku bilgisini, hedef imgenin renk bilgisi ile birleştirerek yeni bir imge oluşturmaktadır.

| $\mathbf{T}$: Hedef İmge |  $\mathbf{S}$: Kaynak İmge | $\Omega$: Kaynak Maskesi | $\mathbf{I}$: Birleştirilmiş İmge
:-------:|:----:|:----:|:---:|:---:|
![Poisson Image Editing][bear_background] | ![Poisson Image Editing][bear_foreground] | ![Poisson Image Editing][bear_foreground_mask] | ![Poisson Image Editing][bear_blended_image] |
![Poisson Image Editing][turkey_background] | ![Poisson Image Editing][turkey_foreground] | ![Poisson Image Editing][turkey_foreground_mask] | ![Poisson Image Editing][turkey_blended_image] |
![Poisson Image Editing][tshirt_background] | ![Poisson Image Editing][tshirt_foreground] | ![Poisson Image Editing][tshirt_foreground_mask] | ![Poisson Image Editing][tshirt_blended_image] |

Yazıda yer alan analizlerin yapıldığı kod parçaları, görseller ve kullanılan veri setlerine [poisson_image_editing](https://github.com/cescript/poisson_image_editing) GitHub sayfası üzerinden erişebilirsiniz.

**Referanslar**
* Pérez, Patrick, Michel Gangnet, and Andrew Blake. "Poisson image editing." ACM SIGGRAPH 2003 Papers. 2003. 313-318.

[RESOURCES]: # (List of the resources used by the blog post)
[bear_background]: resources//bear_background.png
[bear_foreground]: resources//bear_foreground.png
[bear_foreground_mask]: resources//bear_foreground_mask.png
[bear_foreground_mask_explained]: resources//bear_foreground_mask_explained.png
[bear_blended_image]: resources//bear_blended_image.png
[poisson_image_editing_border]: resources//poisson_image_editing.svg

[turkey_background]: resources//turkey_background.png
[turkey_foreground]: resources//turkey_foreground.png
[turkey_foreground_mask]: resources//turkey_foreground_mask.png
[turkey_blended_image]: resources//turkey_blended_image.png

[tshirt_background]: resources//tshirt_background.png
[tshirt_foreground]: resources//tshirt_foreground.png
[tshirt_foreground_mask]: resources//tshirt_foreground_mask.png
[tshirt_blended_image]: resources//tshirt_blended_image.png
